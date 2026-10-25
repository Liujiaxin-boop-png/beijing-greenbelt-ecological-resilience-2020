"""Compare a real raster rerun with the archived numerical study references.

References are used only here, after run_full.py has finished. These checks test
reproduction of the stated computation, not the ecological truth of its proxies.
"""
from pathlib import Path
import json,csv,math,hashlib,argparse
import numpy as np
import rasterio

B=Path(__file__).resolve().parent
def readcsv(path):
    with Path(path).open(encoding='utf-8-sig',newline='')as f:return list(csv.DictReader(f))
def number(value):
    try:return float(value)
    except(TypeError,ValueError):return float('nan')
def js(path,value):Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def verify(output):
    log=json.loads((output/'run_summary.json').read_text(encoding='utf-8'))
    computed=readcsv(output/'tables/grid_components_all_scenarios_scales.csv')
    reference=readcsv(B/'data/reference/archived_components_all_scenarios_scales.csv')
    key=lambda r:(r['grid_id'],r['hq_scenario'])
    cur={key(r):r for r in computed};ref={key(r):r for r in reference}
    if set(cur)!=set(ref):raise ValueError('Grid/scenario identity set differs from archive')
    fields=['R','Re','A','A4','ER_arithmetic','ER_geometric','ER_adjacency4','HQ_known_fraction','Re_known_fraction']
    errors={f:0. for f in fields};missing={f:0 for f in fields};checked=0
    for ident,r in cur.items():
        for f in fields:
            a,b=number(r[f]),number(ref[ident][f]);checked+=1
            if math.isfinite(a)!=math.isfinite(b):missing[f]+=1
            elif math.isfinite(a):errors[f]=max(errors[f],abs(a-b))
    grid_pass=all(v<1e-11 for v in errors.values())and not any(missing.values())
    towns=readcsv(output/'tables/township_means_35_two_products.csv')
    tidx={(r['town_name'],r['scenario']):r for r in towns}
    tref=readcsv(B/'data/reference/township_four_models_reference_35.csv')
    town_errors={};town_rank_errors=[]
    for product,scenario in [('CLCD','CLCD_mid'),('WC','WC_receiver_mid')]:
        for variant in ['arithmetic','geometric']:
            f='ER_'+variant;col=product+'_'+variant+'_area_weighted_mean';rankcol=product+'_'+variant+'_descending_rank_of_35'
            town_errors[col]=max(abs(number(tidx[(r['town_name'],scenario)][f])-number(r[col]))for r in tref)
            vals={r['town_name']:number(tidx[(r['town_name'],scenario)][f])for r in tref}
            for r in tref:
                # Historical town means have no ties; report exactly that ranking.
                rank=1+sum(v>vals[r['town_name']]for v in vals.values())
                if rank!=int(r[rankcol]):town_rank_errors.append({'town':r['town_name'],'model':col,'new':rank,'reference':int(r[rankcol])})
    diagref=readcsv(B/'data/reference/township_diagnostics_reference_35.csv')
    diagfields=['R','Re','A','ER_arithmetic','ER_geometric','ER_CRITIC_arithmetic_diagnostic','ER_arith_omit_R','ER_arith_omit_Re','ER_arith_omit_A','area_km2']
    diagnostic_errors={f:max(abs(number(tidx[(r['town_name'],'CLCD_mid')][f])-number(r[f]))for r in diagref)for f in diagfields}
    quads=readcsv(output/'tables/quadrant_means_two_products.csv');qidx={(r['quadrant'],r['scenario']):r for r in quads}
    qref=readcsv(B/'data/reference/quadrant_four_models_reference.csv');quad_errors={}
    for product,scenario in [('CLCD','CLCD_mid'),('WC','WC_receiver_mid')]:
        for variant in ['arithmetic','geometric']:
            f='ER_'+variant;col=product+'_'+variant+'_area_weighted_mean'
            quad_errors[col]=max(abs(number(qidx[(r['quadrant'],scenario)][f])-number(r[col]))for r in qref)
    raster_checks=[]
    fingerprints=B/'verification/expected_raster_value_fingerprints.json'
    if not fingerprints.exists():raise FileNotFoundError('Frozen original-raster fingerprint ledger is missing')
    for r in json.loads(fingerprints.read_text(encoding='utf-8'))['rasters']:
        path=output/'rasters'/r['recomputed_name']
        with rasterio.open(path)as ds:
            arr=ds.read(1);digest=hashlib.sha256(arr.tobytes(order='C')).hexdigest()
            ok=(digest==r['value_array_sha256'] and list(ds.shape)==r['shape'] and list(ds.transform)[:6]==r['affine'] and ds.crs==rasterio.crs.CRS.from_wkt(r['crs_wkt']) and ds.nodata==r['nodata'])
        raster_checks.append({'file':str(path.relative_to(B)).replace('\\','/'),'same_values_mask_shape_transform_crs_as_archive':bool(ok),'value_array_sha256':digest})
    all_pass=grid_pass and not town_rank_errors and all(v<1e-10 for v in town_errors.values())and all(v<1e-10 for v in diagnostic_errors.values())and all(v<1e-10 for v in quad_errors.values())and all(r['same_values_mask_shape_transform_crs_as_archive']for r in raster_checks)
    result={'status':'PASS_ACTUAL_RASTER_RERUN_MATCHES_ARCHIVED_NUMERICS'if all_pass else'FAIL_OR_DIFFERENCE_REQUIRES_REVIEW','prepared_raster_rerun':True,'scenario_grid_rows_checked':len(computed),'individual_grid_values_checked':checked,'grid_max_absolute_errors':errors,'grid_missingness_disagreements':missing,'township_count':len(tref),'town_four_model_max_absolute_errors':town_errors,'town_four_model_ranking_disagreements':town_rank_errors,'town_main_and_diagnostic_max_absolute_errors':diagnostic_errors,'quadrant_four_model_max_absolute_errors':quad_errors,'raster_count_checked':len(raster_checks),'raster_checks':raster_checks,'recomputed_k':log['computed_k_from_main_D'],'recomputed_Dmax':log['main_reference_Dmax'],'CRITIC_weights':log['CRITIC_weights'],'all_input_files_unchanged':log['all_inputs_sha256_unchanged'],'basis':'Prepared rasters -> six HQ/D scenarios -> full-square R/Re/A and 4-neighbour alternative -> arithmetic/geometric/CRITIC -> exact-area town/quadrant statistics. Historical outputs used only for verification.','not_claimed':['Reacquisition or complete preprocessing of original remote-sensing products','Independent ecological-process validation','Byte-identical GIS/Office/PDF/PNG figure files','Public redistribution permission for all bundled data'],'calculation_output_directory':str(output.relative_to(B)).replace('\\','/')}
    js(B/'verification/full_rerun_verification.json',result)
    print(json.dumps({k:v for k,v in result.items()if k!='raster_checks'},ensure_ascii=False,indent=2))
    if not all_pass:raise SystemExit(2)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='outputs/full_recomputed');a=p.parse_args();verify(B/a.output)
