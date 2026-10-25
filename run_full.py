"""Recompute the 2020 study from the bundled aligned input rasters.

The numerical path is independent of data/reference. It follows the study's
pinned nearest-source HQ equations, full-square aggregation and exact stored
reporting-area weights. It does not run the official InVEST application.
"""
from pathlib import Path
import argparse,csv,json,hashlib,math,time,datetime,sys
import numpy as np
import rasterio
from rasterio.transform import from_origin
from scipy import ndimage
from src.core import exact_pixel_square,cohesion_class,aggregate_components,critic_weights

BASE=Path(__file__).resolve().parent
def readcsv(p):
    with Path(p).open(encoding='utf-8-sig',newline='')as f:return list(csv.DictReader(f))
def savecsv(p,rows):
    with Path(p).open('w',encoding='utf-8-sig',newline='')as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,obj):Path(p).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def log(msg):print(datetime.datetime.now(datetime.timezone.utc).isoformat()+' '+msg,flush=True)
def finite(x):return x is not None and math.isfinite(float(x))

def run(config):
    started=time.perf_counter();cfg=json.loads(config.read_text(encoding='utf-8'))
    output=BASE/cfg['output_dir']
    if output.exists()and any(output.iterdir()):raise FileExistsError('Choose a new empty output_dir; previous runs are not overwritten.')
    for d in [output,output/'rasters',output/'tables']:d.mkdir(parents=True,exist_ok=True)
    input_manifest=json.loads((BASE/'input_manifest.json').read_text(encoding='utf-8'))
    input_checks=[]
    for r in input_manifest['inputs']:
        ok=sha(BASE/r['file'])==r['sha256'];input_checks.append({'file':r['file'],'sha256_matches':ok})
        if not ok:raise ValueError('Input hash mismatch: '+r['file'])
    source_rasters={};nodata={};profile=None
    paths={**cfg['inputs'],'worldcover':cfg['worldcover'],'k_reference_mask':cfg['k_reference_mask'],**{'safe_'+s:p for s,p in cfg['safe_support_by_distance'].items()}}
    for key,path in paths.items():
        with rasterio.open(BASE/path)as ds:
            if profile is None:profile=ds.profile.copy();aff=ds.transform;crs=ds.crs;shape=ds.shape
            if(ds.transform!=aff or ds.crs!=crs or ds.shape!=shape):raise ValueError('Input alignment differs: '+key)
            source_rasters[key]=ds.read(1);nodata[key]=ds.nodata
    assert aff.a==30 and aff.e==-30 and aff.b==aff.d==0
    assert list(aff)[:6]==cfg['expected_raster_affine']
    assert crs==rasterio.crs.CRS.from_wkt((BASE/cfg['analysis_crs_wkt']).read_text(encoding='utf-8'))
    sens=readcsv(BASE/cfg['sensitivity'])
    lut={key:np.full(256,np.nan,dtype='float32')for key in ['habitat','built','road']}
    for row in sens:
        for key in lut:lut[key][int(row['lucode'])]=float(row[key])
    coefficients={int(r['lucode']):float(r['coefficient'])for r in readcsv(BASE/cfg['coefficients'])}
    grids=readcsv(BASE/cfg['grid_bounds'])
    for g in grids:
        for key in ['scale_m','ix','iy']:g[key]=int(g[key])
        for key in ['left','bottom','right','top','grid_area_m2','town_area_m2','belt_area_m2']:g[key]=float(g[key])
    worldcover=np.zeros(shape,'uint8')
    mapping={10:2,20:3,30:4,40:1,50:8,60:7,70:6,80:5,90:9}
    for old,new in mapping.items():worldcover[source_rasters['worldcover']==old]=new
    products={'CLCD':source_rasters['lulc'],'WorldCover':worldcover}
    distances={};threat_positive_counts={}
    for key in ['built','road']:
        a=source_rasters[key];positive=(a!=0)&(a!=nodata[key])
        if not positive.any():raise ValueError('No positive threat source: '+key)
        threat_positive_counts[key]=int(positive.sum())
        log('Computing nearest-positive-source distance: '+key)
        distances[key]=ndimage.distance_transform_edt(~positive).astype('float32')
    log('Computing land-cover composition and 4/8-neighbour cohesion from rasters')
    landscape={}
    for product,land in products.items():
        for g in grids:
            bounds=tuple(g[k]for k in ['left','bottom','right','top']);area=g['grid_area_m2']
            local,wt,rr,cc=exact_pixel_square(land,aff,bounds)
            total=float(wt.sum());amounts={j:float(wt[local==j].sum())for j in range(1,10)}
            known_area=sum(amounts[j]for j in coefficients)
            Re=sum(amounts[j]*v for j,v in coefficients.items())/known_area if known_area and abs(known_area-area)<1e-4 else float('nan')
            centred=local[np.ix_(rr,cc)]
            complete=bool(np.all((centred>0)&(centred<=9)))and abs(total-area)<1e-4
            target=np.isin(centred,cfg['grid']['bluegreen_classes'])
            A=cohesion_class(target,8)[0]if complete else float('nan')
            A4=cohesion_class(target,4)[0]if complete else float('nan')
            landscape[(g['grid_id'],product)]={'Re':Re,'A':A,'A4':A4,'Re_known_fraction':known_area/area}
        log('Landscape metrics complete: '+product)
    allrows=[];rasters=[];fixed_k=None;actual_dmax=None;reference_count=None
    raw_weights=cfg['hq']['raw_weights'];weight_sum=sum(raw_weights.values());z=cfg['hq']['z']
    raster_profile=profile.copy();raster_profile.update(dtype='float32',nodata=-9999,count=1,compress='deflate',predictor=3,tiled=True,blockxsize=256,blockysize=256)
    for scenario,multiplier in cfg['distance_multipliers'].items():
        safe=source_rasters['safe_'+scenario]==1
        for product,land in products.items():
            ident=('CLCD_'if product=='CLCD'else'WC_receiver_')+scenario
            valid=safe&np.isfinite(lut['habitat'][land]);D=np.zeros(shape,'float32')
            for key in ['built','road']:
                cutoff=cfg['hq']['threat_distance_m'][key]*multiplier/abs(aff.a)
                dist=distances[key];decay=np.where(dist>cutoff,0,(cutoff-dist)/cutoff).astype('float32')
                D+=(decay*lut[key][land]*(raw_weights[key]/weight_sum)).astype('float32')
            if ident=='CLCD_mid':
                reference=valid&(source_rasters['k_reference_mask']==1)
                reference_count=int(reference.sum());actual_dmax=float(D[reference].max());fixed_k=actual_dmax/2
                if abs(fixed_k-cfg['hq']['k'])>1e-12:raise ValueError('Recomputed k differs from frozen main reference')
            if fixed_k is None:raise ValueError('Main CLCD_mid must be evaluated first to freeze k')
            Q=(lut['habitat'][land]*(1-D**z/(D**z+fixed_k**z))).astype('float32')
            for metric,a in [('D',D),('HQ',Q)]:
                result=np.where(valid,a,-9999).astype('float32');path=output/'rasters'/f'{metric}_{ident}.tif'
                with rasterio.open(path,'w',**raster_profile)as ds:
                    ds.write(result,1);ds.update_tags(role='Independent study equation implementation; not official InVEST application',scenario=ident,metric=metric,unknown='NoData -9999; not zero',year=2020)
                rasters.append({'scenario':ident,'metric':metric,'file':str(path.relative_to(BASE)).replace('\\','/'),'value_array_sha256':hashlib.sha256(result.tobytes(order='C')).hexdigest(),'valid_cells':int(valid.sum())})
            quality=np.where(valid,Q,-9999).astype('float32')
            for g in grids:
                loc,wt,_,_=exact_pixel_square(quality,aff,tuple(g[k]for k in ['left','bottom','right','top']))
                good=np.isfinite(loc)&(loc!=-9999)&(loc>=0)&(loc<=1)
                support=float(wt[good].sum());area=g['grid_area_m2']
                R=float(np.sum(loc[good]*wt[good])/support)if support and abs(support-area)<1e-4 else float('nan')
                ll=landscape[(g['grid_id'],product)]
                allrows.append({'grid_id':g['grid_id'],'scale_m':g['scale_m'],'ix':g['ix'],'iy':g['iy'],'hq_scenario':ident,'product':product,'town_area_m2':g['town_area_m2'],'belt_area_m2':g['belt_area_m2'],'grid_area_m2':area,'HQ_known_fraction':support/area,'R':R,'Re':ll['Re'],'A':ll['A'],'A4':ll['A4'],'Re_known_fraction':ll['Re_known_fraction']})
            log('HQ rasters and grid integration complete: '+ident)
    primary=[r for r in allrows if r['scale_m']==1000 and r['hq_scenario']=='CLCD_mid'and r['town_area_m2']>0]
    assert len(primary)==2478
    X=np.array([[r[k]for k in ['R','Re','A']]for r in allrows]);Xp=np.array([[r[k]for k in ['R','Re','A']]for r in primary])
    assert np.isfinite(Xp).all();cw=critic_weights(Xp)
    alt4=X.copy();alt4[:,2]=[r['A4']for r in allrows]
    computed={'ER_arithmetic':aggregate_components(X),'ER_geometric':aggregate_components(X,mode='geometric'),'ER_adjacency4':aggregate_components(alt4,mode='geometric'),'ER_arith_adjacency4':aggregate_components(alt4),'ER_CRITIC_arithmetic_diagnostic':aggregate_components(X,cw),'ER_arith_omit_R':aggregate_components(X[:,[1,2]]),'ER_arith_omit_Re':aggregate_components(X[:,[0,2]]),'ER_arith_omit_A':aggregate_components(X[:,[0,1]])}
    for i,r in enumerate(allrows):
        for key,value in computed.items():r[key]=float(value[i])
    savecsv(output/'tables/grid_components_all_scenarios_scales.csv',allrows)
    # Save the real main maps at their native reporting grid resolution.
    for scale in cfg['scales_m']:
        subset=[r for r in allrows if r['scale_m']==scale and r['hq_scenario']=='CLCD_mid']
        nx=max(r['ix']for r in subset)+1;ny=max(r['iy']for r in subset)+1
        for metric in ['R','Re','A','ER_geometric','ER_arithmetic']:
            a=np.full((ny,nx),-9999,np.float32)
            for r in subset:
                if finite(r[metric]):a[ny-1-r['iy'],r['ix']]=r[metric]
            path=output/'rasters'/f'CLCD_mid_{scale}m_{metric}.tif'
            with rasterio.open(path,'w',driver='GTiff',height=ny,width=nx,count=1,dtype='float32',crs=crs,transform=from_origin(cfg['grid']['origin_xy'][0],cfg['grid']['origin_xy'][1]+ny*scale,scale,scale),nodata=-9999,compress='deflate')as ds:
                ds.write(a,1);ds.set_band_description(1,metric)
            rasters.append({'scenario':'CLCD_mid','scale_m':scale,'metric':metric,'file':str(path.relative_to(BASE)).replace('\\','/'),'value_array_sha256':hashlib.sha256(a.tobytes(order='C')).hexdigest()})
    # Full 35-town and four-quadrant values use exact precomputed area weights.
    fields=['R','Re','A','ER_arithmetic','ER_geometric','ER_CRITIC_arithmetic_diagnostic','ER_arith_omit_R','ER_arith_omit_Re','ER_arith_omit_A']
    def aggregate_context(parts,keyfields):
        result=[]
        for scenario in ['CLCD_mid','WC_receiver_mid']:
            index={r['grid_id']:r for r in allrows if r['scale_m']==1000 and r['hq_scenario']==scenario}
            totals={}
            for part in parts:
                r=index.get(part['grid_id'])
                if r is None:continue
                area=float(part['intersection_area_m2']);key=tuple(part[k]for k in keyfields)
                v=totals.setdefault(key,{'area':0.,**{f:0. for f in fields}});v['area']+=area
                for f in fields:v[f]+=area*r[f]
            for key,v in sorted(totals.items()):
                result.append({**dict(zip(keyfields,key)),'scenario':scenario,'area_km2':v['area']/1e6,**{f:v[f]/v['area']for f in fields}})
        return result
    towns=aggregate_context(readcsv(BASE/cfg['intersections']),['town_name','town_code'])
    quadrants=aggregate_context(readcsv(BASE/cfg['quadrant_intersections']),['quadrant'])
    savecsv(output/'tables/township_means_35_two_products.csv',towns)
    savecsv(output/'tables/quadrant_means_two_products.csv',quadrants)
    summaries=[]
    for scale in cfg['scales_m']:
        for scenario in [product+s for s in cfg['distance_multipliers']for product in ['CLCD_','WC_receiver_']]:
            subset=[r for r in allrows if r['scale_m']==scale and r['hq_scenario']==scenario and r['town_area_m2']>0]
            for metric in ['R','Re','A','ER_arithmetic','ER_geometric','ER_arith_adjacency4','ER_CRITIC_arithmetic_diagnostic']:
                good=[r for r in subset if finite(r[metric])];area=sum(r['town_area_m2']for r in good)
                if not area:continue
                mu=sum(r['town_area_m2']*r[metric]for r in good)/area
                summaries.append({'scale_m':scale,'scenario':scenario,'metric':metric,'grid_count':len(good),'area_km2':area/1e6,'area_weighted_mean':mu,'zero_area_percent':sum(r['town_area_m2']for r in good if r[metric]==0)/area*100})
    savecsv(output/'tables/area_weighted_summaries.csv',summaries)
    source_same=all(sha(BASE/r['file'])==r['sha256']for r in input_manifest['inputs'])
    assert source_same
    summary={'status':'ACTUAL_PREPARED_RASTER_TO_ALL_GRID_AND_TOWN_RESULTS_COMPLETED','python':sys.version.split()[0],'numpy':np.__version__,'scipy':__import__('scipy').__version__,'rasterio':rasterio.__version__,'elapsed_seconds':time.perf_counter()-started,'aligned_raster_shape':list(shape),'config':str(config.relative_to(BASE)).replace('\\','/'),'config_sha256':sha(config),'run_script_sha256':sha(__file__),'grid_context_count':len(grids),'main_scope_1km_count':len(primary),'computed_scenario_grid_rows':len(allrows),'HQ_scenarios':6,'township_count_per_product':35,'computed_k_from_main_D':fixed_k,'main_reference_Dmax':actual_dmax,'k_reference_pixel_centers':reference_count,'CRITIC_weights':dict(zip(['R','Re','A'],map(float,cw))),'component_Pearson_unweighted':np.corrcoef(Xp,rowvar=False).tolist(),'positive_threat_cell_counts':threat_positive_counts,'all_inputs_sha256_unchanged':source_same,'input_hash_checks':input_checks,'raster_outputs':rasters,'official_InVEST_application_used':False,'independent_ecological_validation':False,'scope_limits':[cfg['prepared_input_boundary'],'Original boundary polygon generation and raw CLCD/WorldCover/OSM acquisitions are upstream, not rerun here.','No claim that GIS/Excel/Word/PDF figure binary files are reproduced identically.','Boundary-derived context and raster redistribution permissions require author review before uploading data publicly.']}
    js(output/'run_summary.json',summary);log('Completed real prepared-raster calculation; see run_summary.json')
    return summary

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',default='config/full_2020.json');args=parser.parse_args()
    run(BASE/args.config)
