#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 15:26:25 2025

@author: Benedikt Meylahn

With these functions we calculate the exposure in the CBS 500mx500m square data for the paper.
"""
import re
import numpy as np
import math
import geopandas as gpd
import pandas as pd

def distance(center_coord, coord, dim):
    ec, nc = map(int, re.findall(r'-?\d+', center_coord))
    eo, no = map(int, re.findall(r'-?\d+', coord))
    e_dist = abs(ec-eo)#*dim/100 #we are scaling down by 100m for numerical stability
    n_dist = abs(nc-no)#*dim/100 #we are scaling down by 100m for numerical stability
    temp = np.hypot(e_dist, n_dist)
    #print(temp)
    return temp
    
def get_neighbours(coord, data, dim):
    '''
    Old Code for concentric squares, now we use the manhattan distance. kept in case comparison desirable.
    '''
    e_str, n_str = re.findall(r'-?\d+', coord)
    e = int(e_str)
    n = int(n_str)

    e_width = len(e_str)
    n_width = len(n_str)

    neighbours = []
    pointer = int(dim/100)
    for i in [-pointer, 0, pointer]:
        for j in [-pointer, 0, pointer]:
            if i == 0 and j == 0:
                continue
            new_e = str(e+i).zfill(e_width)
            new_n = str(n+j).zfill(n_width)
            new_coord = f"E{new_e}N{new_n}"

            if new_coord in data.index:
                neighbours.append(new_coord)

    return neighbours
    
def get_Man_neighbours(coord, data, dim):
    e_str, n_str = re.findall(r'-?\d+', coord)
    e = int(e_str)
    n = int(n_str)

    e_width = len(e_str)
    n_width = len(n_str)

    neighbours = []
    pointer = int(dim/100)
    for i, j in [(-pointer,0), (pointer,0), (0,-pointer), (0,pointer)]:
        new_e = str(e+i).zfill(e_width)
        new_n = str(n+j).zfill(n_width)
        new_coord = f"E{new_e}N{new_n}"

        if new_coord in data.index:
            neighbours.append(new_coord)

    return neighbours

def clean_percentage(p):
    """
    Returns:
        (value, is_valid)

    value   : float in [0, 1] if valid
    is_valid: False if missing, NaN, or suppressed (negative)
    """

    if p is None:
        return None, False

    try:
        p = float(p)
    except (TypeError, ValueError):
        return None, False

    if math.isnan(p) or p < 0:
        return None, False

    return p / 100.0, True

def local_composition(residents, P_dutch, P_foreigners, P_NL_nonEU, P_NL_EU, P_notNL_EU):
    '''
    Parameters
    ----------
    residents : ?
        Total number of residents in the cell
        data.at[coord, 'aantaal_inwoners'].
    P_dutch : TYPE
        percentage born in the Netherlands, Dutch heritage
        data.at[coord, 'percentage_nederlands'].
    P_foreigners : TYPE
        percentage of residents born outside of Europe
        data.at[coord, 'percentage_niet_nederlands'].
    P_NL_nonEU 
        percentage of residents born in NL but non-EU heritage
        data.at[coord, 'percentage_NL_nonEU'
    P_NL_EU 
        percentage of residents born in NL but non-Dutch, EU heritage
        data.at[coord, 'percentage_NL_EU'  
    P_notNL_EU 
        percentage of residents born outside NL but with EU heritage
        data.at[coord, 'percentage_notNL_EU'
    Returns
    -------
    imputed where necessary number of residents, number of dutch residents and number of foreign residents:
        Res, Dutch, Foreigners.

    '''
    try:
        residents = float(residents)
        residents = max(0,residents)
    except:
        residents = 0
    
    #interrup here if there are zero inhabitants
    if residents == 0:
        return 0, 0, 0, 0, 0
    else:
        
        flt_dutch, valid_dutch = clean_percentage(P_dutch)
        flt_foreign, valid_foreign = clean_percentage(P_foreigners)
        
        flt_NL_nonEU , valid_NL_nonEU = clean_percentage(P_NL_nonEU)
        flt_NL_EU , valid_NL_EU = clean_percentage(P_NL_EU)
        flt_notNL_EU, valid_notNL_EU = clean_percentage(P_notNL_EU)
        
        
        known_shares = []
        for val, valid in [
            (flt_NL_nonEU, valid_NL_nonEU),
            (flt_NL_EU, valid_NL_EU),
            (flt_notNL_EU, valid_notNL_EU)]:
            
            if valid:
                known_shares.append(val)

        resid = sum(known_shares)
        
        if valid_dutch and valid_foreign :
            # check sum bigger than 1 + rounding error. If so, treat dutch as true
            #if flt_dutch + flt_foreign > 1.001:
            #    flt_foreign = 1-flt_dutch
                
            # now we know all values acceptable; proceed assignment without imputation
            num_dutch = residents * flt_dutch
            num_foreign = residents * flt_foreign
            
            return residents, num_dutch, flt_dutch, num_foreign, flt_foreign
        
        elif (not valid_dutch) and valid_foreign:
            # % Dutch suppressed but % foreign, not. Impute that all who are not foreign, are dutch.
            num_foreign = residents * flt_foreign
            num_resid = residents * resid
            num_dutch = residents - num_foreign - num_resid
            flt_dutch = 1-flt_foreign - resid
            
            return residents, num_dutch, flt_dutch, num_foreign, flt_foreign
        
        elif valid_dutch and (not valid_foreign) :
            # % foreign suppressed, but % dutch not. Impute that all who are not dutch, are foreign.
            num_dutch = residents * flt_dutch
            num_resid = residents * resid
            num_foreign = residents - num_dutch - num_resid
            flt_foreign = 1-flt_dutch - resid
            
            return residents, num_dutch, flt_dutch, num_foreign, flt_foreign
        
        else:
            # Neither % NL or % foreign reported, all (non residuals) are Dutch. Immigration is low nationally
            num_dutch = residents - residents*resid
            num_foreign = 0
            flt_dutch = 1.0 - resid
            flt_foreign = 0.0
            
            return residents, num_dutch, flt_dutch, num_foreign, flt_foreign
        
def imm_expos(center_coord, data, max_count, dim, dist_constant, power_factor):
    '''
    Important to note the max_count, is on population individuals reached. Furthermore this is not exactly the same as in Brown and Enos 2021 (https://doi.org/10.1038/s41562-021-01066-z). The aggregated nature of the data requires us to complete a whole layer of cells. Thus the layer in which max_count population is reached is the last layer included.
    '''
    visited = set()
    frontier = [center_coord]
    tot_exp = 0
    tot_weight = 0
    tot_iso = 0
    count = 0

    residents, dutchies, p_dutch, buit, p_buit = local_composition(
       residents =  data.at[center_coord, 'aantal_inwoners'], 
       P_dutch = data.at[center_coord,'percentage_nederlands'], 
       P_foreigners = data.at[center_coord, 'percentage_niet_nederlands'], 
       P_NL_nonEU = data.at[center_coord, 'percentage_NL_nonEU'], 
       P_NL_EU = data.at[center_coord, 'percentage_NL_EU'],
       P_notNL_EU = data.at[center_coord, 'percentage_notNL_EU'])
    

    while frontier != [] and count < max_count:
        next_frontier = []

        if count > max_count:
            break
        for coord in frontier:
            if coord in visited:
                continue
                
            visited.add(coord)
            dist = distance(center_coord, coord, dim)
            # the dist_constant is measured in 100m, due the scaling in the distance calculator, consistent across data
            weight = 1/((dist+dist_constant)**power_factor)
            
            local_res, local_dutch, local_P_dutch, local_foreign, local_P_foreign = local_composition(
                residents =  data.at[coord, 'aantal_inwoners'], 
                P_dutch = data.at[coord,'percentage_nederlands'], 
                P_foreigners = data.at[coord, 'percentage_niet_nederlands'], 
                P_NL_nonEU = data.at[coord, 'percentage_NL_nonEU'], 
                P_NL_EU = data.at[coord, 'percentage_NL_EU'],
                P_notNL_EU = data.at[coord, 'percentage_notNL_EU'])
 
            
            tot_weight += max(0,weight*local_res)
            tot_exp += max(0,weight*local_foreign)
            tot_iso += max(0,weight*local_dutch)
            count += max(0,local_res)
            

            # add neighbours to our next frontier:
            #for neigh in get_neighbours(coord, data, dim): # old concentric squares
            for neigh in get_Man_neighbours(coord, data, dim): # new manhattan diamond
                if neigh not in visited:
                    next_frontier.append(neigh)
        ''' 
        printout to see progress without overloading the console:
        if count >= max_count:
            if np.random.rand()<0.0001:
                print(f"Cell finished. Current count: {count}")'''
        frontier = next_frontier

    if tot_weight == 0:
        return 0, 0, dutchies, p_buit, residents, buit, 0, 0
    else:
        return tot_exp/tot_weight, (tot_exp/tot_weight)*dutchies, dutchies, p_buit, residents, buit, tot_iso/tot_weight, (tot_iso/tot_weight)*dutchies

def municipalities_closeby(geoDF, municipality_name, distance_thresh):
    '''
        This helper is intended to create lists of statcodes (municipality codes)
        which are close to some specified municipality. This is important for 
        fixed effects modelling later on in the analysis.
        The expected input geopandas dataframe is the gemeentengrenzen:
            path_gem_shape_2023 = 'Grid_data/cbsgebiedsindelingen2016-2025/cbsgebiedsindelingen2023.gpkg'
            gemeentegrenzen = gpd.read_file(path_gem_shape_2023, layer = 'gemeente_gegeneraliseerd')
        
        It returns a modified gemeentengrenzen gdf which makes future calcs quicker. The second
        item it returns in the main item of interest. The list of statcodes of municipalities
        within the specified distance from the municipality chosen.
    '''
    geoDF = geoDF.copy()
    
    if not geoDF.crs or not geoDF.crs.is_projected:
        raise ValueError("GeoDataFrame must have a projected CRS")
        
    if 'rep_pt' not in geoDF.columns:
        geoDF['rep_pt'] = geoDF.geometry.representative_point()
        
    sel = geoDF['statnaam'] == municipality_name
    if not sel.any():
        raise ValueError(f"Municipality '{municipality_name}' not found")
        
    # muni_pt = geoDF.loc[sel, 'rep_pt'].iloc[0] # old representative point version
    selected_geom = geoDF.loc[sel, 'geometry'].iloc[0]
    if f"dist_{municipality_name}" not in geoDF.columns:
        geoDF[f"dist_{municipality_name}"] = geoDF.geometry.distance(selected_geom)
        # geoDF[f"dist_{municipality_name}"] = geoDF['rep_pt'].distance(muni_pt) # in combo with representative points.
        
    muni_list = geoDF[geoDF[f"dist_{municipality_name}"]<distance_thresh]['statcode'].to_list()

    return geoDF, muni_list

def compute_exposure_gdf(gdf, exposure_k, dimension, const_c, power_factor_a, suffix = ''):
    gdf = gdf.copy()
    
    indices = gdf.index.to_list()
    results = [imm_expos(idx, gdf, max_count = exposure_k, dim = dimension, dist_constant = const_c, power_factor = power_factor_a) for idx in indices]
    
    exposure_df = pd.DataFrame(results, columns=[f'exposure{suffix}', f'exposure_popNL{suffix}', 
                                                 f'popNL{suffix}', f'p_buit{suffix}', f'pop{suffix}', 
                                                 f'pop_buit{suffix}', f'isol{suffix}', f'isol_popNL{suffix}'])
   
    exposure_df.insert(0, f'crs28992res{dimension}m', indices)
    gdf = gdf.drop(columns=[f'exposure{suffix}', f'exposure_popNL{suffix}', 
                                                 f'popNL{suffix}', f'p_buit{suffix}', f'pop{suffix}', 
                                                 f'pop_buit{suffix}', f'isol{suffix}', f'isol_popNL{suffix}'], errors='ignore')

    calculated_gdf = gdf.merge(exposure_df, left_index=True, right_on=f'crs28992res{dimension}m').set_index(f'crs28992res{dimension}m')
    
    del results, exposure_df
    return calculated_gdf

# Population weighted variables:
def population_weighted(gdf, value_col, pop_col, out_prefix, sentinel=None):
    gdf = gdf.copy()
    if value_col not in gdf.columns:
        gdf[value_col] = 0
    s = gdf[value_col]
    if sentinel is not None:
        s = s.replace(sentinel, np.nan)

    valid = s.notna()

    gdf[f'{out_prefix}_pop'] = gdf[pop_col].where(valid, 0)
    gdf[f'{out_prefix}_dot_pop'] = s.fillna(0) * gdf[f'{out_prefix}_pop']

    return gdf

def aggregate_regions(calculated_gdf, agg_level = 'statcode', suffix=''):
    aggregation_schema = {
        f'exposure{suffix}': 'mean',
        f'exposure_popNL{suffix}': 'sum',
        f'exposure_popNL10': 'sum',
        f'exposure_popNL100': 'sum',
        f'exposure_popNL1000': 'sum',
        f'exposure_popNL5000': 'sum',
        f'exposure_popNL10000': 'sum',
        f'isol_popNL{suffix}':'sum',
        f'popNL{suffix}': 'sum',
        f'p_buit{suffix}': 'mean',
        f'pop{suffix}': 'sum',
        f'pop_buit{suffix}': 'sum',
        'exposure': 'mean',
        'exposure_popNL': 'sum',
        'isol_popNL':'sum',
        'popNL': 'sum',
        'p_buit': 'mean',
        'pop': 'sum',
        'pop_buit': 'sum',
        'stedelijk_filt' : 'mean',
        'aantal_part_huishoudens' : 'sum',
        'num_low' : 'sum',
        'num_high' : 'sum',
        'hh_dot_mean' : 'sum',
        'sted_pop' : 'sum',
        'woz_dot_pop': 'sum',
        'woz_pop':'sum',
        'oad_pop':'sum',
        'aantal_inwoners_15_tot_25_jaar':'sum',
        'aantal_inwoners_25_tot_45_jaar':'sum', 
        'aantal_inwoners_45_tot_65_jaar':'sum', 
        'aantal_inwoners_65_jaar_en_ouder':'sum',
        'aantal_uitkering' : 'sum',
        'oad_dot_pop':'sum'
    }
    #treat missing values:
    try:
        calculated_gdf['stedelijk_filt'] = calculated_gdf['stedelijkheid'].replace(-99997, 5) # missing, thus least urban
    except:
        calculated_gdf['stedelijk_filt'] = 0 # fake value for all squares
        
    calculated_gdf['aantal_inwoners_15_tot_25_jaar'] = calculated_gdf['aantal_inwoners_15_tot_25_jaar'].replace(-99997, 0) # suppressed -> small number -> 0
    calculated_gdf['aantal_inwoners_25_tot_45_jaar'] = calculated_gdf['aantal_inwoners_25_tot_45_jaar'].replace(-99997, 0)
    calculated_gdf['aantal_inwoners_45_tot_65_jaar'] = calculated_gdf['aantal_inwoners_45_tot_65_jaar'].replace(-99997, 0)
    calculated_gdf['aantal_inwoners_65_jaar_en_ouder'] = calculated_gdf['aantal_inwoners_65_jaar_en_ouder'].replace(-99997, 0)
    
    calculated_gdf['aantal_part_huishoudens'] = calculated_gdf['aantal_part_huishoudens'].replace(-99997, 0)
    
    calculated_gdf['aantal_uitkering'] = calculated_gdf['aantal_uitkering'].replace(-99997, np.nan)
    
    calculated_gdf['sted_pop'] = calculated_gdf['stedelijk_filt']*calculated_gdf[f'pop{suffix}']
  
 
    calculated_gdf = population_weighted(
            calculated_gdf,
            value_col='gemiddelde_woz_waarde_woning',
            pop_col=f'pop{suffix}',
            out_prefix='woz', 
            sentinel = [-99997, -99995, -99995.0, '-99995', '99995.0']
            )
   
    calculated_gdf = population_weighted(
            calculated_gdf,
            value_col='omgevingsadressendichtheid',
            pop_col=f'pop{suffix}',
            out_prefix='oad',
            sentinel=[-99997, -99995, -99995.0, '-99995', '99995.0']
            )
    
    valid_agg = {
        col: func
        for col, func in aggregation_schema.items()
        if col in calculated_gdf.columns
    }
    agg_df = calculated_gdf.groupby(agg_level).agg(valid_agg)
    # new per person variables:
    for suff in ['', 10, 100, 1000, 5000, 10000]:
        if f'exposure_popNL{suff}' in agg_df.columns:
            agg_df[f'pers_weigh_exp{suff}'] = agg_df[f'exposure_popNL{suff}'] / agg_df[f'popNL'].replace(0, np.nan)
        if f'isol_popNL{suff}' in agg_df.columns:
            agg_df[f'pers_weigh_isol{suff}'] = agg_df[f'isol_popNL{suff}'] / agg_df[f'popNL'].replace(0, np.nan)
    
    agg_df['voting_pop'] = (agg_df['aantal_inwoners_15_tot_25_jaar'] 
                        +agg_df['aantal_inwoners_25_tot_45_jaar'] 
                         +agg_df['aantal_inwoners_45_tot_65_jaar'] 
                         +agg_df['aantal_inwoners_65_jaar_en_ouder'])
                         
    agg_df['work_pop'] = (agg_df['aantal_inwoners_15_tot_25_jaar'] 
                        +agg_df['aantal_inwoners_25_tot_45_jaar'] 
                         +agg_df['aantal_inwoners_45_tot_65_jaar'])  
                         
    valid = (
        agg_df['aantal_uitkering'].notna() &
        agg_df['work_pop'].notna() &
        (agg_df['work_pop'] > 0)
            )               
    agg_df['share_WW'] = np.nan
    agg_df.loc[valid, 'share_WW'] = (
        agg_df.loc[valid, 'aantal_uitkering'] /
        agg_df.loc[valid, 'work_pop']
            )

    agg_df['share_15_to_25'] = agg_df['aantal_inwoners_15_tot_25_jaar']/agg_df['voting_pop']
    agg_df['share_25_to_45'] = agg_df['aantal_inwoners_25_tot_45_jaar']/agg_df['voting_pop']
    agg_df['share_45_to_65'] = agg_df['aantal_inwoners_45_tot_65_jaar']/agg_df['voting_pop']

    #agg_df['pers_low'] = agg_df['num_low'] / agg_df['aantal_part_huishoudens'].replace(0, np.nan) # (for later once this data arrives)
    #agg_df['pers_high'] = agg_df['num_high'] / agg_df['aantal_part_huishoudens'].replace(0, np.nan)
    #agg_df['mean_income'] = agg_df['hh_dot_mean'] / agg_df['aantal_part_huishoudens'].replace(0, np.nan)

    agg_df['sted_pop_mean'] = agg_df['sted_pop']/agg_df[f'pop{suffix}'].replace(0, np.nan)
 
    agg_df['woz_pop_mean'] = agg_df['woz_dot_pop']/agg_df['woz_pop'].replace(0, np.nan)

    agg_df['oad_pop_mean'] = agg_df['oad_dot_pop']/agg_df['oad_pop'].replace(0, np.nan)

    '''
    agg_df = agg_df[['pers_weigh_exp_new','pers_weigh_isol', f'pop{suffix}', 
                        f'pop_buit{suffix}', f'popNL{suffix}', 'mean_income','pers_high', 
                        'pers_low', 'stedelijk_filt', "sted_pop_mean", 'mean_woz',
                       'share_15_to_25', 'share_25_to_45', 'share_45_to_65', 'oad_mean']]
    '''
    
    agg_df = agg_df.reset_index()
    return agg_df
    
def renaming_keeping(path_to_gdf):
    '''
        Reads a geopandas dataframe into existance and restricts it to the appropriate columns
    '''
    gdf = gpd.read_file(path_to_gdf)
    
    cols = ['crs28992res500m','crs28992res100m', 'geometry', 'aantal_inwoners', 'aantal_part_huishoudens',
                  'percentage_geb_nederland_herkomst_nederland','percentage_geb_buiten_nederland_herkmst_buiten_europa',
                  'percentage_geb_nederland_herkomst_overig_europa', 'percentage_geb_nederland_herkomst_buiten_europa', 
                  'percentage_geb_buiten_nederland_herkomst_europa',
                  'stedelijkheid', 'omgevingsadressendichtheid', 'aantal_personen_met_uitkering_onder_aowlft',
                  'gemiddeld_inkomen_huishouden', 'percentage_laag_inkomen_huishouden','percentage_hoog_inkomen_huishouden', 
                 'aantal_inwoners_15_tot_25_jaar', 'aantal_inwoners_25_tot_45_jaar', 'aantal_inwoners_45_tot_65_jaar', 'aantal_inwoners_65_jaar_en_ouder', 
                 'aantal_woningen', 'gemiddelde_woz_waarde_woning', 'percentage_koopwoningen', 'percentage_huurwoningen']
    cols_to_keep = [c for c in cols if c in gdf.columns]
    
    rename_dict = {'percentage_geb_nederland_herkomst_nederland': 'percentage_nederlands', 
        'percentage_geb_nederland_herkomst_buiten_europa' : 'percentage_NL_nonEU', 
        'percentage_geb_nederland_herkomst_overig_europa' : 'percentage_NL_EU',
        'percentage_geb_buiten_nederland_herkomst_europa': 'percentage_notNL_EU',
            'percentage_geb_buiten_nederland_herkmst_buiten_europa': 'percentage_niet_nederlands', 'aantal_personen_met_uitkering_onder_aowlft' : 'aantal_uitkering'}
    
    if 'crs28992res100m' in cols_to_keep:
        index_name = 'crs28992res100m'
    elif 'crs28992res500m' in cols_to_keep:
        index_name = 'crs28992res500m'
    else:
        raise ValueError('Geopandas dataframe gdf must have either crs28992res500m, or crs28992res100m')
    
    gdf = (
        gdf[cols_to_keep]
        .rename(columns=rename_dict)
        .set_index(index_name)
        )
    return gdf
    
def combine_region_grid(grid_gdf, region_gdf):
    assert grid_gdf.crs.is_projected

    grid_gdf = grid_gdf.copy()
    grid_gdf["pt"] = grid_gdf.geometry.representative_point()
    
    assigned_df = gpd.sjoin(
        grid_gdf.set_geometry("pt"),
        region_gdf,
        how="inner",
        predicate="within"
    ).set_geometry("geometry")
    return assigned_df
    
    
## Election data management, initially only summing the right wing votes

def merging_election_regions(elect_gdf, region_gdf):
    elect_gdf = elect_gdf.drop(columns=['index_right'], errors='ignore')
    region_gdf = region_gdf.drop(columns=['index_right'], errors='ignore')
    
    gdf_joined = gpd.sjoin(
        elect_gdf,
        region_gdf,
        how="left",
        predicate="intersects"   # points within polygons
        )

    pary_cols = ['PVV (Partij voor de Vrijheid)', 'Forum voor Democratie', 'GROENLINKS / Partij van de Arbeid (PvdA)', 'SP (Socialistische Partij)', 
             'JA21', 'Staatkundig Gereformeerde Partij (SGP)' , 'D66', 'DENK', 'BIJ1', 'Partij voor de Dieren' , 'ChristenUnie',
            'Volt', 'VVD', 'CDA','Nieuw Sociaal Contract', 'geldige stembiljetten']

    gdf_joined[pary_cols] = (
        gdf_joined[pary_cols]
          .apply(pd.to_numeric, errors="coerce")
          .fillna(0)
          .astype("Int64")
        )

    elect_agg = gdf_joined.groupby('statcode').agg({'PVV (Partij voor de Vrijheid)':'sum', 'Forum voor Democratie':'sum', 
                                            'GROENLINKS / Partij van de Arbeid (PvdA)':'sum',
                                            'SP (Socialistische Partij)':'sum', 'JA21':'sum', 'Staatkundig Gereformeerde Partij (SGP)' : 'sum',
                                            'D66' : 'sum', 'DENK' : 'sum', 'BIJ1': 'sum', 'Partij voor de Dieren' : 'sum', 'ChristenUnie' : 'sum',
                                                'Volt' : 'sum', 'VVD' : 'sum', 'CDA' : 'sum', 'Nieuw Sociaal Contract' : 'sum', 
                                            'geldige stembiljetten':'sum'})


    elect_agg['RW'] = (elect_agg['PVV (Partij voor de Vrijheid)'] + elect_agg['Forum voor Democratie']+elect_agg['JA21'])/ elect_agg['geldige stembiljetten'].replace(0, np.nan)
    elect_agg['PVV'] = elect_agg['PVV (Partij voor de Vrijheid)']/elect_agg['geldige stembiljetten'].replace(0, np.nan)
    elect_agg = elect_agg.drop(columns = ['PVV (Partij voor de Vrijheid)', 'Forum voor Democratie', 'GROENLINKS / Partij van de Arbeid (PvdA)', 'SP (Socialistische Partij)', 
             'JA21', 'Staatkundig Gereformeerde Partij (SGP)' , 'D66', 'DENK', 'BIJ1', 'Partij voor de Dieren' , 'ChristenUnie',
            'Volt', 'VVD', 'CDA','Nieuw Sociaal Contract'])
    return elect_agg
    
def log_scaling(gdf, col):
    '''
    assumes col is strictly positive!
    '''
    gdf = gdf.copy()
    gdf[f'{col}_log'] = np.log(gdf[col])
    gdf[f'{col}_log_scaled'] = (gdf[f'{col}_log'] -  gdf[f'{col}_log'].mean())/ gdf[f'{col}_log'].std()

    return gdf
    
def zscore_scaling(gdf, col):
    gdf = gdf.copy()
    gdf[f'{col}_z_scaled'] = (gdf[col] -  gdf[col].mean())/ gdf[col].std()

    return gdf
    
def pop_group(x, small_lim=1000):
    if x < small_lim:
        return "small"
    elif x >= 10_000:
        return "large"
    else:
        return "medium"
