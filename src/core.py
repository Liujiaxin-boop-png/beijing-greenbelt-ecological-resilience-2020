"""Transparent ecological-condition proxy calculations; original sources in REFERENCES.md.

Prepared-input implementation of the study equations, not the official InVEST engine.
All paths and acquisition responsibilities belong to the calling configuration.
"""
import math
import numpy as np
from scipy import ndimage

def cohesion_class(binary, adjacency=8):
    """FRAGSTATS class COHESION /100, full-square native pixel centres.
    Empty-target=0 is an explicit extension; missing support must be masked by caller.
    Patches use 8 or 4 adjacency, exposed perimeter always uses 4 cell faces.
    """
    Z=binary.size
    if Z<=1:return float('nan'),0
    structure=np.ones((3,3),int) if adjacency==8 else ndimage.generate_binary_structure(2,1)
    labels,n=ndimage.label(binary,structure=structure)
    if n==0:return 0.,0
    areas=np.bincount(labels.ravel(),minlength=n+1).astype(float)
    per=np.zeros(n+1,float)
    for q in [labels[0,:],labels[-1,:],labels[:,0],labels[:,-1]]:per+=np.bincount(q,minlength=n+1)
    for first,second in [(labels[:-1,:],labels[1:,:]),(labels[:,:-1],labels[:,1:])]:
        different=first!=second
        per+=np.bincount(first[different],minlength=n+1);per+=np.bincount(second[different],minlength=n+1)
    value=(1-per[1:].sum()/(per[1:]*np.sqrt(areas[1:])).sum())/(1-1/math.sqrt(Z))
    if value < -1e-10 or value>1+1e-10:raise ValueError('COHESION outside theoretical bounds')
    return float(np.clip(value,0,1)),int(n)

def exact_pixel_square(a,aff,bounds):
    """Pixel/full-square intersection areas; supports a north-up native 30m grid."""
    left,bottom,right,top=bounds
    c0=max(0,int(math.floor((left-aff.c)/aff.a)));c1=min(a.shape[1],int(math.ceil((right-aff.c)/aff.a)))
    r0=max(0,int(math.floor((aff.f-top)/(-aff.e))));r1=min(a.shape[0],int(math.ceil((aff.f-bottom)/(-aff.e))))
    cols=np.arange(c0,c1);rows=np.arange(r0,r1);xs=aff.c+cols*aff.a;ys=aff.f+rows*aff.e
    wx=np.maximum(0,np.minimum(xs+aff.a,right)-np.maximum(xs,left))
    wy=np.maximum(0,np.minimum(ys,top)-np.maximum(ys+aff.e,bottom));wt=wy[:,None]*wx[None,:]
    return a[r0:r1,c0:c1],wt,(ys+aff.e/2<=top)&(ys+aff.e/2>bottom),(xs+aff.a/2>=left)&(xs+aff.a/2<right)

def habitat_quality(lulc,threat_arrays,threat_nodata,safe,sensitivity,distances_m,raw_weights,pixel_size=30,k=.39925000071525574,z=2.5):
    """Nearest-positive-threat EDT + linear attenuation, float32 like archived run.
    Distances are nearest-source distances, not a convolution of every threat cell.
    """
    lut={key:np.full(256,np.nan,dtype='float32') for key in ['habitat','built','road']}
    for row in sensitivity:
        for key in lut:lut[key][int(row['lucode'])]=float(row[key])
    valid=(safe==1)&np.isfinite(lut['habitat'][lulc]);D=np.zeros(lulc.shape,'float32');weights=sum(raw_weights.values())
    for key in ['built','road']:
        a=threat_arrays[key];nd=threat_nodata[key];positive=(a!=0)&(a!=nd)
        if not positive.any():raise ValueError('No positive valid threat source: '+key)
        distance=ndimage.distance_transform_edt(~positive).astype('float32');cutoff=distances_m[key]/pixel_size
        decay=np.where(distance>cutoff,0,(cutoff-distance)/cutoff).astype('float32')
        D+=(decay*lut[key][lulc]*(raw_weights[key]/weights)).astype('float32')
    Q=(lut['habitat'][lulc]*(1-D**z/(D**z+k**z))).astype('float32')
    return D,Q,valid

def aggregate_components(values,weights=None,mode='arithmetic'):
    a=np.asarray(values,float)
    if weights is None:weights=np.ones(a.shape[-1])/a.shape[-1]
    weights=np.asarray(weights,float)
    if np.any(weights<0) or not np.isclose(weights.sum(),1):raise ValueError('Weights must be nonnegative and sum to one')
    if mode=='arithmetic':return a@weights
    if mode!='geometric':raise ValueError(mode)
    positive=weights>0;out=np.exp(np.log(np.maximum(a[:,positive],np.finfo(float).tiny))@weights[positive])
    out[np.any((a==0)&(weights>0),axis=1)]=0;out[np.any(~np.isfinite(a)&(weights>0),axis=1)]=np.nan
    return out

def critic_weights(values):
    """Diagnostic on existing native 0-1 component scales, not minmax-standardized.
    This data-dispersion weight is not evidence of ecological importance.
    """
    a=np.asarray(values,float);a=a[np.isfinite(a).all(axis=1)]
    corr=np.corrcoef(a,rowvar=False);sigma=a.std(axis=0,ddof=1);information=sigma*np.sum(1-corr,axis=1)
    if not np.isfinite(information).all() or information.sum()<=0:raise ValueError('Undefined CRITIC weights')
    return information/information.sum()

def weighted_townships(rows,intersections,fields):
    index={r['grid_id']:r for r in rows};totals={}
    for part in intersections:
        rec=index.get(part['grid_id'])
        if rec is None:continue
        key=(part['town_name'],part.get('town_code',''));weight=float(part['intersection_area_m2'])
        if weight<=0:continue
        out=totals.setdefault(key,{'area_m2':0.,**{f:0. for f in fields}});out['area_m2']+=weight
        for field in fields:out[field]+=weight*float(rec[field])
    return [dict(town_name=key[0],town_code=key[1],area_km2=out['area_m2']/1e6,**{f:out[f]/out['area_m2'] for f in fields}) for key,out in sorted(totals.items())]
