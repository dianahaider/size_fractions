#!/usr/bin/env python

#for importing, formatting and data manipulation
import pandas as pd
import numpy as np
import glob
import tempfile
from qiime2 import Artifact
import yaml
import json
import os
import re

#for plotting
import matplotlib, random
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from matplotlib.patches import Patch
import seaborn as sns
from collections import Counter
from collections import defaultdict

#sns.set(style="whitegrid")
import plotly.express as px
from IPython.display import display
from IPython.display import HTML
from sklearn.ensemble import IsolationForest

from pandas.plotting import register_matplotlib_converters
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
register_matplotlib_converters()
import scipy as sp
import statsmodels.api as sm

#for statistical analyses
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from skbio.stats.distance import permanova
from skbio import DistanceMatrix
from scipy.spatial.distance import cdist
from skbio.stats.composition import clr
from skbio.stats.composition import ancom
from skbio.diversity.alpha import shannon
import scipy.stats as stats
import statsmodels.api as sa
import statsmodels.formula.api as sfa
import scikit_posthocs as spPH
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score, ConfusionMatrixDisplay
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from scipy.stats import randint

from plyer import notification

# Special thanks to Alex Manuele https://github.com/alexmanuele
def consolidate_tables(MG, frac=None):
    if MG == '16S' or  MG == 'chloroplast':
        comm = '02-PROKs'
    elif MG == '18S' :
        comm = '02-EUKs'
    print('Community is '+comm)

    if frac==None:
        table_list = glob.glob('{0}/dada2_SF/table.qza'.format('/Users/Diana/Documents/escuela/phd/size_fractions/BB22_size-fraction-comparison-analysed/to_transfer/'+comm))
        print("Found all "+MG+" tables.")
    else:
        table_list = glob.glob('{0}/*/table.qza'.format('/Users/Diana/Documents/escuela/phd/size_fractions/BB22_size-fraction-comparison-analysed/to_transfer/'+comm))
        print("Found all "+MG+" tables.")


    dataframes = []
    for table_path in table_list:
        with tempfile.TemporaryDirectory() as tempdir:
            #load table, dump contents to tempdir
            table = Artifact.load(table_path)
            #Make sure the tables are all FeatureFrequency type
            assert str(table.type) == 'FeatureTable[Frequency]', "{0}: Expected FeatureTable[Frequency], got {1}".format(table_path, table.type)
            Artifact.extract(table_path, tempdir)
            #get the provenance form the tempdir and format it for DF
            prov = '{0}/{1}/provenance/'.format(tempdir, table.uuid)
            action = yaml.load(open("{0}action/action.yaml".format(prov), 'r'), Loader=yaml.BaseLoader)
            paramlist = action['action']['parameters']
            paramlist.append({'table_uuid': "{}".format(table.uuid)})
            paramdict = {}
            for record in paramlist:
                paramdict.update(record)

            # Get the data into a dataframe
              #Biom data
            df = table.view(pd.DataFrame).unstack().reset_index()
            df.columns = ['feature_id', 'sample_name', 'feature_frequency']
            df['table_uuid'] = ["{}".format(table.uuid)] * df.shape[0]
              #param data
            pdf = pd.DataFrame.from_records([paramdict])
              #merge params into main df
            df = df.merge(pdf, on='table_uuid')


            #I like having these columns as the last three. Makes it more readable
            cols = df.columns.tolist()
            reorder = ['sample_name', 'feature_id', 'feature_frequency']
            for val in reorder:
                cols.append(cols.pop(cols.index(val)))
            df = df[cols]
            df['table_path'] = [table_path] * df.shape[0]
            df['sample_name'] = df['sample_name'].str.replace('-', '.')
            dataframes.append(df)

            # Adding table_id, forward and reverse trim columns
            #df['table_id'] = str(table_path.split('/')[-3]) #add a table_id column
            #df['forward_trim'], df['reverse_trim'] = df['table_id'].str.split('R', 1).str
            #df['forward_trim'] = df['forward_trim'].map(lambda x: x.lstrip('F'))
            #df["forward_trim"] = pd.to_numeric(df["forward_trim"])
            #df["reverse_trim"] = pd.to_numeric(df["reverse_trim"])

    #Stick all the dataframes together
    #outputfile="merged_all_tables.tsv"
    df = pd.concat(dataframes)
    df['sample_name'] = df['sample_name'].str.replace(r'\.S([1-9]|[1-9][0-9]|[1-9][0-9][0-9]).L001\.','', regex=True)
    df['sample_name'] = df['sample_name'].str.replace('V4V5.','', regex=True)
    df['sample_name'] = df['sample_name'].str.replace(r'\.LandS.S([1-9]|[1-9][0-9]|[1-9][0-9][0-9]).L001','P', regex=True)
    #df.to_csv(comm+'/merged_all_tables.tsv', sep='\t', index=False)
    print("Successfully saved all tables.")

    if MG == 'chloroplast':
        comm = 'chloroplast'

    return df, comm



def merge_metadata(df, all_md):
    #df = pd.read_csv('02-PROKs/'+'/merged_all_tables.tsv', sep='\t')

    depth_num = {
        "A": 1,
        "B": 5,
        "C": 10,
        "D": 60,
        "E": 30
    }

    tables = df[['sample_name', 'feature_id', 'feature_frequency']].copy()
    tables.rename(columns={'sample_name':'sampleid'}, inplace=True)

    all_md['sampleid'] = all_md['sampleid'].str.replace('_', '.')
    merged = pd.merge(tables,all_md, on='sampleid', how='left') #all_md is the metadata file
    merged = merged[merged.feature_frequency != 0]

    merged['year'] = 2022

    merged["size_code"] = merged["sampleid"].str.extract(r'[1-9][0-9]?[A-E]([L-S])')
    merged["size_code"] = merged["size_code"].fillna('W')
    merged["depth_code"] = merged["sampleid"].str.extract(r'[1-9][0-9]?([A-E])')
    merged['depth']= merged['depth_code'].map(depth_num)
    merged["weekn"] = merged["sampleid"].str.extract(r'\.([1-9][0-9]?)[A-E]')
    merged['weekn'] = pd.to_numeric(merged['weekn'])
    merged['depth'] = pd.to_numeric(merged['depth'])
    merged['date'] = merged.groupby('weekn', as_index=False)['date'].transform('first')

    merged['Total'] = merged['feature_frequency'].groupby(merged['sampleid']).transform('sum')
    merged['ratio'] = merged['feature_frequency']/merged['Total']
    merged['nASVs'] = merged['feature_id'].groupby(merged['sampleid']).transform('count')
    merged['weekdepth'] = merged["weekn"].astype(str) + merged["depth"].astype(str)
    merged['avg'] = merged['nASVs'].groupby(merged['weekdepth']).transform('mean')
    merged['diff'] = merged['nASVs'] - merged['avg']

    print('Set up metadata ...')

    #merged.to_csv(comm+'/merged_asvs_metadata.tsv', sep = '\t')
    print('Saved merged_asvs_metadata.tsv')

    return merged

def make_defract(all_md, separated):
    #make sure all size codes are indicated
    all_md["size_code"] = all_md["sampleid"].str.extract(r'[1-9][0-9]?[A-E]([L-S])')
    all_md["size_code"] = all_md["size_code"].fillna('W')

    #only keep values from weeks 1 to 16
    sep_SL = all_md[
        (all_md.size_code != "W") &
        (all_md.size_code != "P")
    ]
    sep_SL = sep_SL.drop(sep_SL[sep_SL.weekn > 16].index)

    #sum [DNA] of small and large size fractions
    sep_SL['[DNAt]'] = sep_SL.groupby(['weekn', 'depth'])['[DNA]ng/ul'].transform('sum')

    #separate small and large size fraction
    sep_S = sep_SL[sep_SL.size_code == 'S']
    sep_L = sep_SL[sep_SL.size_code == 'L']

    #calculate DNA proportion per size fraction
    sep_SL['DNApr'] = sep_SL['[DNA]ng/ul']/sep_SL['[DNAt]']

    #merge with separated on common columns to get corresponding rel. abundances
    sep_SL = sep_SL[['sampleid', 'DNApr', '[DNAt]']].copy()
    sepSLRA = pd.merge(separated, sep_SL, on=['sampleid'], how='left') #all_md is the metadata file

    #exclude ASVs from the whole water
    sep_SLRA = sepSLRA[
        (sepSLRA.size_code != 'W') &
        (sepSLRA.size_code != 'P')
    ]

    #calculate corrected per sample ratio, and corrected feature frequency of de-fractionated samples
    sep_SLRA['Newfeature_frequency'] = sep_SLRA['feature_frequency'] * sep_SLRA['DNApr']
    sep_SLRA['Newff'] = sep_SLRA.groupby(['feature_id', 'weekn', 'depth'])['Newfeature_frequency'].transform('sum')


    #sep_SLRA = sep_SLRA.drop(['sampleid', 'size_code'], axis=1)
    sep_SLRA['sampleid'] = "BB22." + sep_SLRA['weekn'].astype(str) + sep_SLRA['depth_code'] + "SL"

    #uncomment the line below if keeping small and large original sample
    #sep_SLRA['size_code'] = sep_SLRA['size_code'] + '-DFr'

    #uncomment the line above if merging smallandlarge
    sep_SLRA['size_code'] = 'SL'

    #drop unecessary columns which might rise merging conflicts
    sep_SLRA = sep_SLRA.drop(['feature_frequency', 'Total', 'ratio', 'nASVs', 'weekdepth', 'avg',
                              'diff', 'extraction_date', '[DNA]ng/ul', 'A260/280', 'A260/230',
                              'Newfeature_frequency'], axis=1)
    sep_SLRA.rename(columns={'Newff':'feature_frequency'}, inplace=True)
    sep_SLRA = sep_SLRA.drop_duplicates()

    #recalculate ratios
    sep_SLRA['Total'] = sep_SLRA['feature_frequency'].groupby(sep_SLRA['sampleid']).transform('sum')
    sep_SLRA['ratio'] = sep_SLRA['feature_frequency']/sep_SLRA['Total']
    sep_SLRA['nASVs'] = sep_SLRA['feature_id'].groupby(sep_SLRA['sampleid']).transform('nunique')

    sep_SLRA = sep_SLRA.drop_duplicates()

    #make new df dependingg on plotting needs
    sep_WO = separated[separated.size_code == "W"]
    sep_WO = sep_WO.drop_duplicates()

    sep_PO = separated[separated.size_code == "P"]
    sep_PO = sep_PO.drop_duplicates()

    sep_S = separated[separated.size_code == "S"]
    sep_L = separated[separated.size_code == "L"]


    sep_WO.reset_index(inplace=True, drop=True)
    sep_SLRA.reset_index(inplace=True, drop=True)

    #newseparated = pd.concat([sep_SLRA.reset_index(drop=True), sep_WO.reset_index(drop=True)], axis=0).reset_index(drop=True)
    newseparated = pd.concat([sep_SLRA, sep_WO, sep_PO, sep_L, sep_S], ignore_index=True)

    newseparated['weekdepth'] = newseparated["weekn"].astype(str) + newseparated["depth"].astype(str)
    newseparated['avg'] = newseparated['nASVs'].groupby(newseparated['weekdepth']).transform('mean')
    newseparated['diff'] = newseparated['nASVs'] - newseparated['avg']

    newseparated["rank"] = newseparated.groupby("sampleid")["ratio"].rank(method="average", ascending=False)
    newseparated["ranktot"] = newseparated['rank'] / newseparated['nASVs']

    #calculate shannon diversity index
    grouped = newseparated.groupby('sampleid')['feature_frequency'].apply(list)
    diversity = grouped.apply(shannon)
    newseparated['shannon_diversity'] = newseparated['sampleid'].map(diversity)

    #calculate dnaconc of SL (debugged by chatGPT)
    cleaned_data = newseparated.drop_duplicates(subset=["weekn", "depth", "date", "size_code", "[DNA]ng/ul"])

    sl_fill_values = (
        cleaned_data[cleaned_data["size_code"].isin(["S", "L"])]
        .groupby(["weekn", "depth", "date"], as_index=False)["[DNA]ng/ul"]
        .sum()
    )
    sl_fill_values["size_code"] = "SL"
    newseparated.loc[
        (newseparated["size_code"] == "SL") & (newseparated["[DNA]ng/ul"].isna()),
        "[DNA]ng/ul",] = newseparated.merge(
        sl_fill_values,
        on=["weekn", "depth", "date", "size_code"],
        how="left")["[DNA]ng/ul_y"]

    return newseparated

def outlier(df):
    new_df = df.copy()
    numeric_cols = ['L', 'S', 'SL', 'W']

    q1 = np.percentile(new_df[numeric_cols],25, axis=0)
    q3 = np.percentile(new_df[numeric_cols],75, axis=0)
    IQR = q3 - q1
    lower_limit = q1 - (1.5*IQR)
    upper_limit = q3 + (1.5*IQR)
    mask = (new_df[numeric_cols] < lower_limit) | (new_df[numeric_cols] > upper_limit)
    new_df[numeric_cols] = new_df[numeric_cols].mask(mask)
    return new_df


def pick_metadata(comm, merged, depth='all', size_fraction='both', year='all', R='all', F='all', txsubset = 'all'):
#make df of features/composition+run+comm

    depth = depth
    year = year
    size_fraction = size_fraction
    txsubset = txsubset

    files = glob.glob('{0}/*/class/*/data/taxonomy.tsv'.format('/Users/Diana/Documents/escuela/phd/size_fractions/BB22_size-fraction-comparison-analysed/to_transfer'))
    taxos = []
#    if not os.path.exists(path+composition):
#        os.mkdir(path+composition)
    for filename in files:
        tax = pd.read_csv(filename, sep='\t')
        taxos.append(tax)

    print('Appended all taxonomies to taxos')
    taxos = pd.concat(taxos)
    taxos = taxos.rename(columns={"Feature ID": "feature_id"}, errors="raise")
    taxos = taxos.drop_duplicates()

    separated = merged.merge(taxos, how='left', on='feature_id') #merged excludes features of frequency = 0
    separated = separated.drop_duplicates()

    if depth != 'all':
        separated = separated[separated["depth"] == depth]
    if size_fraction != 'both':
        separated = separated[separated["size_fraction"] == size_fraction]

    separated[['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']] = separated['Taxon'].str.split('; ', expand=True)
    cols = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    for col in cols:
        separated[col] = separated[col].fillna('Unassigned')

    separated['Month'] = separated['date'].str.split('-').str[1]

    #separated['total'] = separated.groupby(['table_id','sample-id'])['feature_frequency'].transform('sum')
    #separated['ratio'] = separated['feature_frequency']/(separated['total'])
    #separated_taxonomies = separated.copy()

    #make a dictionary with keys for id-ing the taxon belonging to this sub-community
    #separated_dic = pd.Series(separated.Taxon.values,separated.feature_id.values).to_dict()
    print('Saved separated by metadata dataframe.')
    print('Community is '+comm)

    contaminants = None
    if comm in ['02-PROKs', '02-EUKS']:
        searchfor = ["Cyanobacteria", "Chloroplast"]
        contaminants = separated[separated.Taxon.str.contains('|'.join(searchfor))]
        separated = separated[~separated.Taxon.str.contains('|'.join(searchfor))]
        separated = separated.reset_index(drop=True)
        print('Removed cyanobacteria and chloroplast from '+comm)

        #re-calculate ratios when removing chloroplast for 16S or 18S
        separated['Total'] = separated['feature_frequency'].groupby(separated['sampleid']).transform('sum')
        separated['ratio'] = separated['feature_frequency']/separated['Total']
        separated['nASVs'] = separated['feature_id'].groupby(separated['sampleid']).transform('count')
        separated['weekdepth'] = separated["weekn"].astype(str) + separated["depth"].astype(str)
        separated['avg'] = separated['nASVs'].groupby(separated['weekdepth']).transform('mean')
        separated['diff'] = separated['nASVs'] - separated['avg']

    elif comm == 'chloroplast':
        #run these lines to switch to chloroplast comm
        searchfor = ["Cyanobacteria", "Chloroplast"]
        contaminants = separated[separated.Taxon.str.contains('|'.join(searchfor))]
        separated = contaminants.copy()
        separated = separated.reset_index(drop=True)
        comm = 'chloroplast'
        print('Switched to cyanobacteria and chloroplast')

        #add phytorep taxonomy
        cp_tax = pd.read_csv('chloroplast/taxonomy.tsv', sep='\t')
        cp_tax = cp_tax.rename(columns={"Feature ID": "feature_id", "Taxon": "PRTaxon", "Confidence":"PRConfidence"})
        separated = pd.merge(separated, cp_tax, on="feature_id", how="left")
        separated['PRSpecies'] = separated['PRTaxon'].str.split('|').str[-1]

        #re-calculate ratios when removing chloroplast for 16S
        separated['Total'] = separated['feature_frequency'].groupby(separated['sampleid']).transform('sum')
        separated['ratio'] = separated['feature_frequency']/separated['Total']
        separated['nASVs'] = separated['feature_id'].groupby(separated['sampleid']).transform('count')
        separated['weekdepth'] = separated["weekn"].astype(str) + separated["depth"].astype(str)
        separated['avg'] = separated['nASVs'].groupby(separated['weekdepth']).transform('mean')
        separated['diff'] = separated['nASVs'] - separated['avg']

    return separated, contaminants


def SRA_pairs(comm, SFX, SFY, separated, outliers='None', view=False):

    depths = [1,5,10,30,60]

    df_results = pd.DataFrame(columns = ['Depth', 'X', 'Y', 'Coeff', 'Pvalue', 'Rsq'],
                             index = depths)
#outliers has to be a list of indices to remove from separated
    if outliers != 'None':
        cleaned = separated.drop(outliers, axis=0)
    else:
        cleaned = separated.copy()

    for depth in depths:
        d1 = cleaned.loc[cleaned['depth'] == depth]
        forpl = d1[['ratio', 'feature_id', 'sampleid', 'weekn', 'depth', 'size_code', 'Phylum', 'Family']].copy()
        slwplot = forpl.pivot_table(index=["feature_id", "depth", 'weekn','Phylum', 'Family'], columns="size_code", values='ratio').fillna(0)
        slwplot = slwplot.reset_index()

        #build model
        Y = slwplot[SFY]
        X = slwplot[SFX]

        X = sm.add_constant(X)

        #fit the model
        model = sm.OLS(Y, X, missing='drop')
        model_result = model.fit()
        model_result.summary()

        sns.histplot(model_result.resid);
        mu, std = stats.norm.fit(model_result.resid)
        print('Residuals For depth ' + str(depth) + ': mu=' + str(mu) + 'std=' + str(std))

        if view == True:
            fig, ax = plt.subplots()
            # plot the residuals
            sns.histplot(x=model_result.resid, ax=ax, stat="density", linewidth=0, kde=True)
            ax.set(title="Distribution of residuals", xlabel="residual")

            # plot corresponding normal curve
            xmin, xmax = plt.xlim() # the maximum x values from the histogram above
            x = np.linspace(xmin, xmax, 100) # generate some x values
            p = stats.norm.pdf(x, mu, std) # calculate the y values for the normal curve
            sns.lineplot(x=x, y=p, color="orange", ax=ax)
        if view == True:
            plt.show()
            sns.boxplot(x=model_result.resid, showmeans=True);
            sm.qqplot(model_result.resid, line='s');

        fig, ax = plt.subplots()
        fig.set_size_inches(6, 5)

        sns.set_style('white')
        fig = sm.graphics.plot_fit(model_result,1, vlines=False, ax=ax)
        ax.set_ylabel("Defractionated")
        ax.set_xlabel("Whole")
        ax.set_title("Fitted values linear regression")

        plt.savefig('outputs/'+comm+'/asv_'+str(depth)+'_RL_'+SFX+SFY+'.png', dpi=200, bbox_inches="tight")
        if view == True:
            plt.show()

        Y_max = Y.max()
        Y_min = Y.min()

        plt.figure(figsize=(5, 4))
        ax = sns.scatterplot(x=model_result.fittedvalues, y=Y)
        ax.set(ylim=(Y_min, Y_max))
        ax.set(xlim=(Y_min, Y_max))
        ax.set_xlabel("Predicted value")
        ax.set_ylabel("Observed value")

        X_ref = Y_ref = np.linspace(Y_min, Y_max, 100)
        plt.plot(X_ref, Y_ref, color='red', linewidth=1)

        plt.savefig('outputs/'+comm+'/asv_'+str(depth)+'_RL_'+SFX+SFY+'.png', dpi=200, bbox_inches="tight")
        if view == True:
            plt.show()

        df_results.loc[depth] = [depth, SFX, SFY, model_result.params[0], model_result.pvalues[0], model_result.rsquared]
        df_results.to_csv('outputs/'+comm+'/RL_results'+SFX+SFY+'.csv', index=False)

    return df_results

def timeseries_fid(comm, newseparated, f_id, scl, depth):
    #if all size codes
    newseparated['weekfid'] = newseparated["weekn"].astype(str) + newseparated["feature_id"].astype(str)
    d_spc = newseparated[newseparated.depth == depth]

    #if only SL and W
    #d_spc = newsep2[newsep2.depth == depth]

    d_spc_c = d_spc[d_spc.feature_id == f_id]
    #d_spc_c = d_spc[d_spc.Genus == scl]

    sizecodes = ['S', 'L', 'W', 'SL', 'P']
    palette_colors = sns.color_palette()
    palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}

    sel_cols = d_spc_c[['weekn', 'ratio', 'size_code', 'weekfid']].copy()
    sel_cols.drop_duplicates(inplace=True)
    sel_cols["weekn"] = pd.to_numeric(sel_cols["weekn"])
    #make sure each size code has 16weeks
    sel_cols = (sel_cols.set_index(["size_code", "weekn"])
                .reindex(pd.MultiIndex.from_product([sel_cols["size_code"].unique(), [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]],
                                                    names=["size_code", "weekn"]))
                .reset_index()[sel_cols.columns])
    #make sure each week has all fractions
    sel_cols = (sel_cols.set_index(["weekn","size_code"])
                .reindex(pd.MultiIndex.from_product([sel_cols["weekn"].unique(), ["S", "L", "SL", "W"]],
                                                    names=["weekn","size_code"]))
                .reset_index()[sel_cols.columns])
    #fill with 0
    sel_cols.fillna(0, inplace=True)


    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)

    #g = sns.lmplot(x="weekn", y="ratio", data=sel_cols, truncate=True,
    #           lowess=True, ci=None, scatter_kws={"s": 80}, hue='size_code', palette=palette_dict)


    sel_cols['W_per_wk'] = sel_cols['ratio'].groupby(sel_cols['weekfid']).transform('mean')
    sel_cols['diff_p_id'] = sel_cols['ratio'] - sel_cols['W_per_wk']

    g = sns.lineplot(x="weekn", y="ratio", data=sel_cols, hue='size_code',
                     palette=palette_dict, marker='o')

    if g.get_legend() is not None:
        sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
    else:
        print("No legend attached to the plot.")

    g.set_xticks(np.arange(1, 16.5, 1))
    g.figure.set_size_inches(5,3.5)

    if "__" in scl:
        plt.title(scl.split('__')[1] +', '+str(depth)+ 'm')
    else:
        plt.title(scl +', '+str(depth)+ 'm')

    # Set x-axis label
    plt.xlabel('Time (weeks)')
    # Set y-axis label
    plt.ylabel('Relative abundance')
    #plt.ylabel('Rank')

    plt.savefig('outputs/'+comm+'/D'+str(depth)+scl+f_id+'_lineplot.png', dpi=200, bbox_inches="tight")

def taxbarplot(comm, table, level, depth, topn, colrow): #separated is the df, #level is a string of taxonomic level column name, depth is an integer
    if comm == 'chloroplast':
        level = 'PRTaxon'

    table.loc[table[level].isin(['Unassigned', 'g__uncultured', 's__uncultured']), level] = table['feature_id']


    #get a list of top taxa to provide the palette for the visualisation
    toptaxa = table[['feature_frequency', 'Taxon', 'size_code', 'depth','weekn', level]].copy()
    toptaxa = toptaxa.drop_duplicates()
    df_agg = toptaxa.groupby(['size_code',level, 'depth']).agg({'feature_frequency':sum})
    topd = df_agg['feature_frequency'].groupby(['size_code', 'depth'], group_keys=False).nlargest(topn)
    topd = topd.to_frame()
    topd = topd.reset_index()
    listoftop = topd[level].unique()

    #set a palette for the toptaxa
    hex_colors_dic = {}
    rgb_colors_dic = {}
    hex_colors_only = []
    for name, hex in matplotlib.colors.cnames.items():
        hex_colors_only.append(hex)
        hex_colors_dic[name] = hex
        rgb_colors_dic[name] = matplotlib.colors.to_rgb(hex)

    palette_dict = {taxon: color for taxon, color in zip(listoftop, px.colors.sequential.Plasma)}

    sfd=table[table.depth==depth]
    toptaxa = sfd[['feature_frequency', 'Taxon', 'size_code', 'depth','weekn', level]].copy()
    toptaxa = toptaxa.drop_duplicates()
    df_agg = toptaxa.groupby(['size_code',level, 'depth']).agg({'feature_frequency':sum})
    topd = df_agg['feature_frequency'].groupby('size_code', group_keys=False).nlargest(topn)
    topd = topd.to_frame()
    topd = topd.reset_index()


    df_agg = df_agg.reset_index()
    df_agg['set_name'] = df_agg['size_code']+df_agg['depth'].astype(str)

    cumulab = table[['feature_frequency', 'depth', 'size_code', level]].copy()
    cumulab1 = cumulab.groupby([level]).agg({'feature_frequency':sum})

    resultpivot = df_agg.pivot_table(index=level, columns='set_name', values='feature_frequency')
    resultpivot = resultpivot.fillna(0)
    resultpivot[resultpivot != 0] = 1
    tosave = pd.merge(resultpivot, cumulab1, left_index=True, right_index=True)
    tosave.to_csv('outputs/'+comm+'/relab'+level+'_'+str(depth)+'.csv')

    top10d_list = topd[level].unique()
    top10d = sfd.copy()
    top10d.loc[~top10d[level].isin(top10d_list), level] = 'Other' #isnot in top list
    phyld = top10d.groupby(['size_code','weekn', level])['ratio'].sum()
    phyld = phyld.reset_index()


    fig = px.bar(phyld, x="size_code", y="ratio", facet_col="weekn", color=level, labels={
                     "feature_frequency": "Relative abundance",
                     "size_code": "",
                     "weekn": "w"}, color_discrete_map=palette_dict,
    hover_data=[level])
    fig.update_xaxes(type='category', dtick=1)
    fig.update_layout(
        #title= text="Relative abundance of top"+str(topn) + level + 'observed at Depth' + str(depth),
        bargap=0.1,
        yaxis_title="Relative abundance",
        xaxis_title="Size fraction",
        legend_title=level,
        font=dict(size=8),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )

    fig.show()
    #fig.write_image("outputs/fig1.png")
    #fig.to_image(format="png")

    return phyld, top10d


# In[ ]:


def pcaplot(separated, depth, comm, columnperm, spc, colrow):

    if depth == 'all':
        df = separated.copy()
    else:
        df=separated[separated.depth==depth]


    if 'SL' in separated['size_code'].unique():
        #sizecode palette codes
        sizecodes = ['S', 'L', 'W', 'SL', 'P']
        palette_colors = sns.color_palette()
        palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}
        dicsc = pd.Series(df.size_code.values,index=df.sampleid).to_dict()
        color_rows_sc = {k: palette_dict[v] for k, v in dicsc.items()}
        seriescr = pd.Series(color_rows_sc)

    else:
        #sizecode palette codes
        sizecodes = ['S', 'L', 'W', 'P']
        palette_colors = sns.color_palette()
        palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}
        dicsc = pd.Series(df.size_code.values,index=df.sampleid).to_dict()
        color_rows_sc = {k: palette_dict[v] for k, v in dicsc.items()}
        seriescr = pd.Series(color_rows_sc)

    #month palette code
    if colrow == 'Month':
        df['Month'] = df['date'].str.split('-').str[1]
        months = ['Jan', 'Feb', 'Mar', 'May', 'Apr']
        palette_colors = sns.color_palette("flare")
        palette_dict_month = {monthname: color for monthname, color in zip(months, palette_colors)}
        dic = pd.Series(df.Month.values,index=df.sampleid).to_dict()
        color_rows_month = {k: palette_dict_month[v] for k, v in dic.items()}
        seriesmonthcr = pd.Series(color_rows_month)
        dfcolors = pd.DataFrame({'Month': seriesmonthcr,'Size code':seriescr})

    else:
        df['weekn2'] = df['weekn'].astype(str)
        weeks = list(df['weekn'].unique())
        weeks.sort()
        weeks = [str(x) for x in weeks]
        palette_colors = sns.color_palette("flare", 16)
        palette_dict_weekn = {weekname: color for weekname, color in zip(weeks, palette_colors)}
        dic = pd.Series(df.weekn2.values,index=df.sampleid).to_dict()
        color_rows_weekn = {k: palette_dict_weekn[v] for k, v in dic.items()}
        seriesweekn = pd.Series(color_rows_weekn)
        dfcolors = pd.DataFrame({'Weekn': seriesweekn, 'Size code': seriescr})


    topiv = df[['feature_id', 'feature_frequency', 'sampleid']].copy()
    topiv = topiv.drop_duplicates()

    sfdpiv= topiv.pivot(index='sampleid', columns='feature_id', values='feature_frequency')
    sfdpiv=sfdpiv.fillna(0)
    sfdclr=sfdpiv.mask(sfdpiv==0).fillna(0.1)
    clr_transformed_array = clr(sfdclr)
    samples = sfdpiv.index
    asvs = sfdpiv.columns

    #Creating the dataframe with the clr transformed data, and assigning the sample names
    clr_transformed = pd.DataFrame(clr_transformed_array, columns=asvs)
    #Assigning the asv names
    clr_transformed['samples'] = samples
    clr_transformed = clr_transformed.set_index('samples')
    clr_transformed.head()

    #calculate distance matrix
    dist = cdist(clr_transformed, clr_transformed, 'euclid')
    distance_matrix = pd.DataFrame(dist, columns=samples)
    distance_matrix['samples'] = samples
    distance_matrix = distance_matrix.set_index('samples')

    #format for pca
    dm = DistanceMatrix(distance_matrix)

    pca = PCA(n_components=2)
    pca_features = pca.fit_transform(distance_matrix)

    ####
    sns.set(rc={"figure.figsize":(4, 3)})
    sns.set_style("whitegrid", {'axes.grid' : False})
    plot_df = pd.DataFrame(data = pca_features, columns = ['dim1', 'dim2'], index = sfdpiv.index)
    plot_df['dim1'] = plot_df['dim1']/1000
    plot_df['dim2'] = plot_df['dim2']/1000
    if depth =='all':
        plot_df2 = pd.merge(plot_df,df[['sampleid','size_code','depth']],on='sampleid', how='left')
    else:
        plot_df2 = pd.merge(plot_df,df[['sampleid','size_code','weekn']],on='sampleid', how='left')


    ##divide into pre-post bloom
    def get_stage(weekNb):
        if weekNb < 8:
            return 'Pre-bloom'
        elif weekNb >= 8:
            return 'Bloom'

    if depth != 'all':
        plot_df2['Time'] = plot_df2['weekn'].apply(get_stage)

    plot_df2 = plot_df2.rename(columns={'size_code': 'Size code'})

    pc1v = round(pca.explained_variance_ratio_[0]*100)
    pc2v = round(pca.explained_variance_ratio_[1]*100)

    #plot_df2 = plot_df2.drop_duplicates()
    #dfperm = plot_df2.set_index('sampleid')

    #permanova2 = permanova(dm, dfperm, columnperm)
    #results = permanova2(999)

    #plot

    if depth == 'all':
        var2 = 'depth'
    else:
        var2 = 'Time'

    sns.set_style("white")
    ax=sns.scatterplot(x = 'dim1', y = 'dim2', hue= 'Size code', style=var2, data = plot_df2,
                       palette=palette_dict) #, size = 'Week_Group')#,palette=sns.color_palette("dark:salmon_r", as_cmap=True))
    plt.ylabel('PCo2') #(' + str(pc2v) + '% variance explained)')
    plt.xlabel('PCo1') #(' + str(pc1v) +'% variance explained)')
    ax.set_title('Depth ' + str(depth) + 'm', loc='left', weight='bold')
    plt.legend(frameon=False)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    sns.despine()
    plt.savefig('outputs/'+comm+'/D'+str(depth)+spc+'_PCAplot.png', dpi=200, bbox_inches="tight")
    plt.clf()
    plt.cla()
    plt.close()

    print ( "Components = ", pca.n_components_ , ";\nTotal explained variance = ",
      round(pca.explained_variance_ratio_.sum(),5)  )

    print ("Components 1 and 2 are", pca.explained_variance_ratio_)

    # Retrieve Loadings
    loadings = pca.components_

    # Summarize Loadings by Metadata Category
    metadata_groups = plot_df2[var2].unique()
    metadata_contributions = {}

    for group in metadata_groups:
        group_variables = plot_df2.loc[plot_df2[var2] == group, 'sampleid']
        group_loadings = np.abs(loadings[:, [list(distance_matrix.columns).index(var) for var in group_variables]]).mean(axis=1)
        metadata_contributions[group] = group_loadings

    # Visual Representation
    #for group, contributions in metadata_contributions.items():
    #    plt.barh(contributions, group) #range(1, len(contributions) + 1),

    #plt.ylabel('Principal Component')
    #plt.xlabel('Average Loading Contribution')
    #sns.despine()
    #plt.legend(frameon=False)
    #plt.savefig('outputs/'+comm+'/D'+str(depth)+spc+'_PCAplot_brplot.png', dpi=200, bbox_inches="tight")
    #plt.clf()
    #plt.cla()
    #plt.close()


    ##clustermap
    ax = sns.clustermap(distance_matrix, method="complete", cmap='RdBu', annot=True,
               yticklabels=True, row_colors = dfcolors,
               annot_kws={"size": 7}, figsize=(15,12));
    if colrow == 'Month':
        handles1 = [Patch(facecolor=palette_dict_month[key]) for key in palette_dict_month]
        plt.legend(handles1, palette_dict_month, title='Month',
                   bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper left')

    else:
        handles1 = [Patch(facecolor=palette_dict_weekn[key]) for key in palette_dict_weekn]
        plt.legend(handles1, palette_dict_weekn, title='Week',
                bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper left')


    plt.savefig('outputs/'+comm+'/D'+str(depth)+spc+'_clustermap.png', dpi=200, bbox_inches="tight")
    plt.clf()
    plt.cla()
    plt.close()


    plot_df2['PCo 1'] = pc1v
    plot_df2['PCo 2'] = pc2v
    plot_df2.rename(columns={'Size code': 'size_code'}, inplace=True)
    plot_df2.to_csv('R_results/R_testing_vis/'+str(depth)+'_'+comm+'for_R.csv')

    return pca, pca_features, sfdclr, distance_matrix

def dnacon(newseparated, depth='all', includeSL=True):
    if depth != 'all':
        df=newseparated[newseparated.depth==depth]

    sorted_data = newseparated[["weekn", "[DNA]ng/ul", "size_code",'depth',
                                'date']].sort_values(by=["weekn", "size_code", 'depth'])
    if includeSL == False:
        filtered_df = sorted_data[sorted_data["size_code"] != "SL"]
    else:
        filtered_df = sorted_data
    filtered_df.drop_duplicates().sort_values('[DNA]ng/ul')

    sizecodes = ['Small', 'Large', 'Whole', 'Small + Large', 'P']
    palette_colors = sns.color_palette()
    palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}

    filtered_df["size_code"] = filtered_df["size_code"].map({'S': 'Small', 'L': 'Large', 'W': 'Whole', 'SL':'Small + Large'})

    plt.figure(figsize=(12, 6))
    sns.barplot(
        data=filtered_df,
        x="date",
        y="[DNA]ng/ul",
        hue="size_code",
        palette=palette_dict,
        hue_order = ["Whole", "Small + Large", "Small", "Large"]
    )

    plt.xlabel("Date", fontsize=14, weight='bold')
    plt.ylabel("DNA Concentration (ng/µl)", fontsize=14, weight='bold')
    plt.legend(title="Size fraction", fontsize=10, title_fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    # Display the plot
    plt.tight_layout()
    plt.savefig("outputs/grouped_bar_dnaconc"+str(depth)+".png", dpi=300, bbox_inches="tight")
    plt.show()





def rarefy_curve(comm, newseparated):
    if 'SL' in newseparated['size_code'].unique():
        #sizecode palette codes
        sizecodes = ['S', 'L', 'W', 'SL', 'P']
        palette_colors = sns.color_palette()
        palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}
        dicsc = pd.Series(newseparated.size_code.values,index=newseparated.sampleid).to_dict()
        color_rows_sc = {k: palette_dict[v] for k, v in dicsc.items()}
        seriescr = pd.Series(color_rows_sc)

    else:
        #sizecode palette codes
        sizecodes = ['S', 'L', 'W', 'P']
        palette_colors = sns.color_palette()
        palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}
        dicsc = pd.Series(newseparated.size_code.values,index=newseparated.sampleid).to_dict()
        color_rows_sc = {k: palette_dict[v] for k, v in dicsc.items()}
        seriescr = pd.Series(color_rows_sc)

    ax=sns.lineplot(data=newseparated, x="Total", y="nASVs", hue="size_code", marker="o", linestyle=(0, (1, 10)),
                   palette=palette_dict)
    ax.set(xlabel='Sequences per sample', ylabel='Observed richness')

    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    sns.despine()
    plt.clf()
    plt.savefig('outputs/'+comm+'_rarefaction_curve.png', dpi=200, bbox_inches="tight")


    ax=sns.scatterplot(data=newseparated, x="Total", y="nASVs", hue="size_code", marker="o",
                   palette=palette_dict)
    ax.set(xlabel='Sequences per sample', ylabel='Observed richness')

    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    sns.despine()

    plt.savefig('outputs/'+comm+'_rarefaction_curve_noline.png', dpi=200, bbox_inches="tight")


def roll_avg(comm, table, depth, col, rollingavg=4):
    new2 = table.set_index('weekn')
    new2.sort_index(ascending=True, inplace=True)
    new2 = new2[new2['depth'] == depth]
    analysis = new2[[col, 'size_code']].copy()

    sizecodes = ['S', 'L', 'SL', 'W']
    for sizecode in sizecodes:
        df1 = analysis[analysis['size_code'] == sizecode]
        df1.drop('size_code', axis=1, inplace=True)
        df1.drop_duplicates(inplace=True)

        rolling_mean = df1.rolling(rollingavg).mean()
        rolling_std = df1.rolling(rollingavg).std()

        plt.clf()
        plt.plot(df1, color="blue",label="nASVs")
        plt.plot(rolling_mean, color="red", label="Rolling Mean")
        plt.plot(rolling_std, color="black", label = "Rolling Standard Deviation")

        plt.ylabel('nASVs')
        plt.xlabel('Time (week)')
        plt.title('Rolling mean nASVs at depth' + str(depth) + 'm and size' + sizecode , loc='left', weight='bold')
        #plt.legend(frameon=False)

        plt.savefig('outputs/'+comm+'/D'+str(depth)+sizecode+'_rollavg.png', dpi=200, bbox_inches="tight")

    #get one plot of the average of all size fractions
    df1 = analysis.drop('size_code', axis=1) #drop the SF identifier
    df1['mean_col'] = df1.groupby(df1.index)['nASVs'].mean() #take the average of all fractions per week

    df1.drop('nASVs', axis=1, inplace=True)
    df1.drop_duplicates(inplace=True)

    rolling_mean = df1.rolling(rollingavg).mean()
    rolling_std = df1.rolling(rollingavg).std()

    plt.clf()
    plt.plot(df1, color="blue",label="nASVs")
    plt.plot(rolling_mean, color="red", label="Rolling Mean")
    plt.plot(rolling_std, color="black", label = "Rolling Standard Deviation")

    plt.ylabel('nASVs')
    plt.xlabel('Time (week)')
    plt.title('Rolling mean nASVs at depth' + str(depth) + 'm and all fractions', loc='left', weight='bold')
    #plt.legend(frameon=False)

    plt.savefig('outputs/'+comm+'/D'+str(depth)+'_rollavg.png', dpi=200, bbox_inches="tight")

    analysis['weekn'] = analysis.index
    analysis.drop_duplicates(inplace=True)
    plt.clf()
    sizecodes = ['S', 'L', 'W', 'SL']
    palette_colors = sns.color_palette()
    palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}
    ax = sns.lmplot(x ='weekn', y ='nASVs', data = analysis, hue ='size_code', palette=palette_dict, legend=False)
    plt.ylabel('nASVs')
    plt.xlabel('Time (week)')
    plt.title('Trend of nASVs at depth' + str(depth), loc='left', weight='bold')

    plt.savefig('outputs/'+comm+'/D'+str(depth)+'_nasvtrend.png', dpi=200, bbox_inches="tight")



    return pca, pca_features, sfdclr, distance_matrix

def roll_avg(comm, table, depth, col, rollingavg=4):
    new2 = table.set_index('weekn')
    new2.sort_index(ascending=True, inplace=True)
    new2 = new2[new2['depth'] == depth]
    analysis = new2[[col, 'size_code']].copy()

    sizecodes = ['S', 'L', 'SL', 'W']
    for sizecode in sizecodes:
        df1 = analysis[analysis['size_code'] == sizecode]
        df1.drop('size_code', axis=1, inplace=True)
        df1.drop_duplicates(inplace=True)

        rolling_mean = df1.rolling(rollingavg).mean()
        rolling_std = df1.rolling(rollingavg).std()

        plt.clf()
        plt.plot(df1, color="blue",label="nASVs")
        plt.plot(rolling_mean, color="red", label="Rolling Mean")
        plt.plot(rolling_std, color="black", label = "Rolling Standard Deviation")

        plt.ylabel('nASVs')
        plt.xlabel('Time (week)')
        plt.title('Rolling mean nASVs at depth' + str(depth) + 'm and size' + sizecode , loc='left', weight='bold')
        #plt.legend(frameon=False)

        plt.savefig('outputs/'+comm+'/D'+str(depth)+sizecode+'_rollavg.png', dpi=200, bbox_inches="tight")

    #get one plot of the average of all size fractions
    df1 = analysis.drop('size_code', axis=1) #drop the SF identifier
    df1['mean_col'] = df1.groupby(df1.index)['nASVs'].mean() #take the average of all fractions per week

    df1.drop('nASVs', axis=1, inplace=True)
    df1.drop_duplicates(inplace=True)

    rolling_mean = df1.rolling(rollingavg).mean()
    rolling_std = df1.rolling(rollingavg).std()

    plt.clf()
    plt.plot(df1, color="blue",label="nASVs")
    plt.plot(rolling_mean, color="red", label="Rolling Mean")
    plt.plot(rolling_std, color="black", label = "Rolling Standard Deviation")

    plt.ylabel('nASVs')
    plt.xlabel('Time (week)')
    plt.title('Rolling mean nASVs at depth' + str(depth) + 'm and all fractions', loc='left', weight='bold')
    #plt.legend(frameon=False)

    plt.savefig('outputs/'+comm+'/D'+str(depth)+'_rollavg.png', dpi=200, bbox_inches="tight")

    analysis['weekn'] = analysis.index
    analysis.drop_duplicates(inplace=True)
    plt.clf()
    sizecodes = ['S', 'L', 'W', 'SL']
    palette_colors = sns.color_palette()
    palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}
    ax = sns.lmplot(x ='weekn', y ='nASVs', data = analysis, hue ='size_code', palette=palette_dict, legend=False)
    plt.ylabel('nASVs')
    plt.xlabel('Time (week)')
    plt.title('Trend of nASVs at depth' + str(depth), loc='left', weight='bold')

    plt.savefig('outputs/'+comm+'/D'+str(depth)+'_nasvtrend.png', dpi=200, bbox_inches="tight")

#replacement dictionary for unidentified genus ASV for top 20 ASVs
def apply_replacement(df, column_to_map, column_to_update):
    #incl. 18s, 16s, chloro
    replacement_dict = {'d2d1baa1b5487ea6da52cbe74835d92e': 'Unknown Dinophyceae*', #18s
                      '206758a86c52c6b169e5aea9769ee8a3': 'Unknown Dinophyceae*',
                      'f748ff93828e750fa8212b0c15ebb312': 'Unknown Dinophyceae*',
                      'ea8b5924f9a6a3cf089fb3f5c2ce00b1': 'Unknown Dinophyceae*',
                      'ec0cef5115e1a5e72ba718c19f44ac1e': 'Unknown Dinophyceae*',
                      '4ad35c9c4bc721fe9b6a63b0351e2527': 'Unknown Oligotracheae*',
                      '6071d3037004f54c1009530f02f29072': 'Uncultured synidiales*',
                      'd98e045a8124335c3f7422745eb649df': 'Unknown Dinophyceae*',
                      '579e5ddf97287a76bbbc58b1c42ef2f2': 'Unknown Dinophyceae*',
                      '25679bd7ae54946d9d7348b7fde04db4': 'Unknown Dinophyceae*',
                      '43a118280860f972f1ee6c8813ee31ce': 'Uncultured Spirotrichea*',
                      '2949e15fb8d243ac547148c767e37505': 'Uncultured Choreotrichia*',
                      '8ec8ff981ab1011dda5943044363126d': 'Cyclotella*',
                      '334feb949bba92d67ce2d15e5517f30d': 'Pseudo-nitzschia*',
                      '352da949fea628d994cfa85942d0fd95': 'Pseudopedinella*',
                      'a6403cc4635bf5216e10d8724763fb72': 'Guinardia-like*',
                      '6c987094eb76bff568a4499383aa85e6': 'Pseudopyrenomonadales*',
                      '1f6540f210aa24caa4354afe31cf8e7c': 'Guinardia-like*',
                      'd540a1d35a894bc8e3ce9e2413c61dd1': 'Unclassified vampirovibrio*',
                      '7338e9edfef034ec21924ff03d3f4417': 'Synechococcus_CC9902*',
                      '9dbe1d31f92323dcbdcc77903c6bba89': 'Unidentified Ochromonas*',
                      'fcc8273752120d2dc6503688d9986f2e': 'Asteraceae*',
                      'd64d4fd57ebe3985f0fbc3f68fd18aa0': 'Thermosynechococcus*',
                      'f75a7ff2740c7af21f310955d1fe0528' :'Uncultured Thioglobaceae*',
                      '15949af7fbc59962388dce15963f9cec': 'Uncultured Nitrincolaceae*',
                      '5a94578dd1d7cdd039a52f1c7079f874': 'Uncultured Flavobacteriaceae*',
                      'c76526dd2767e5090c0b1e42096ddb92': 'Uncultured Gimesiaceae*',
                      'f98a2cdb11cbbef735b0705b57171c79': 'Uncultured Thecofilosea*',
                      'f3aa3ab8b0d2ae94859675d59169af75': 'Unidentified Pucciniaceae*'}
    # Create a new column based on the mapping
    updated_column = df[column_to_map].map(replacement_dict).fillna(df[column_to_update])

    # Check if changes were made
    if (df[column_to_update] == updated_column).all():
        print("Nothing to change")
    else:
        print("Values were updated")

    # Update the DataFrame column
    df[column_to_update] = updated_column
    return df


def detect_anomalies(metadata, df, dpt, yr=all, month=all):

    sfd=df[df.depth==dpt]

    md_col = sfd[['event_id', metadata, "year", "month"]].copy()
    md_col = md_col[md_col[metadata].notna()]
    if yr != all:
        #mdcol_yr = md_col[md_col.Year == yr]
        mdcol_yr = md_col[md_col['year'].isin(yr)]
    else:
        mdcol_yr = md_col

    if month != all:
        #mdcol_yr = mdcol_yr[mdcol_yr.Month == month]
        mdcol_yr = mdcol_yr[mdcol_yr['month'].isin(month)]

    mdcol_yr = mdcol_yr.drop(columns=['year', "month"])
    mdcol_yr = mdcol_yr.set_index(['event_id'])

    #modelling time
    outliers_fraction = float(.01)
    scaler = StandardScaler()
    np_scaled = scaler.fit_transform(mdcol_yr.values.reshape(-1, 1))
    data = pd.DataFrame(np_scaled)
    # train isolation forest
    model =  IsolationForest(contamination=outliers_fraction)
    model.fit(data)

    #predict data
    mdcol_yr['anomaly'] = model.predict(data)

    # visualization
    fig, ax = plt.subplots(figsize=(10,6))
    a = mdcol_yr.loc[mdcol_yr['anomaly'] == -1, [metadata]] #anomaly
    ax.plot(mdcol_yr.index, mdcol_yr[metadata], color='black', label = 'Normal')
    ax.scatter(a.index,a[metadata], color='red', label = 'Anomaly')
    #plt.axvline(36, ls='--')
    plt.legend()
    plt.show();
    #add axes names

def save_all4_plots(comm,depth,fid,newseparated):
    sizecodes = ['S', 'L', 'W', 'SL']
    palette_colors = sns.color_palette()
    palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}

    sfd=newseparated[newseparated.depth==depth]
    sfd['weekfid'] = sfd["weekn"].astype(str) + sfd["feature_id"].astype(str)
    sfd['avg_p_id'] = sfd['ratio'].groupby(sfd['weekfid']).transform('mean')
    sfd['diff_p_id'] = sfd['ratio'] - sfd['avg_p_id']

    sfd_f=sfd[sfd.feature_id==fid]

    tax = sfd.loc[sfd["feature_id"] == fid].iloc[0]["Taxon"]
    rk = sfd.loc[sfd["feature_id"] == fid].iloc[0]["rank"]

    ttl = sfd_f['Taxon'].iloc[0]

    sns.set(rc={"figure.figsize":(9, 4)})
    sns.set_style("ticks")
    plt.subplot(221)
    ax=sns.lineplot(data=sfd_f, x="weekn", y="diff_p_id", hue="size_code", palette=palette_dict)#, hue="size_code")
    plt.legend(title='Size code')
    ax.legend([],[], frameon=False)

    plt.ylabel('ΔRA per week')
    plt.xlabel('Week number')


    sfd_f['SCfid'] = sfd_f["size_code"].astype(str) + sfd_f["feature_id"].astype(str)
    sfd_f['avg_p_sc'] = sfd_f['ratio'].groupby(sfd_f['SCfid']).transform('mean')
    sfd_f['diff_p_sc'] = sfd_f['ratio'] - sfd_f['avg_p_sc']


    sns.set_style("ticks")
    plt.subplot(222)
    ax=sns.lineplot(data=sfd_f, x="weekn", y="diff_p_sc", hue="size_code", palette=palette_dict)#, hue="size_code")
    plt.legend(title='Size code')
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), frameon=False)

    plt.ylabel('ΔRA per size fraction')
    plt.xlabel('Week number')


    sns.set_style("ticks")
    plt.subplot(223)
    ax = sns.lineplot(data = sfd_f, x='weekn', y='ratio', hue='size_code', palette=palette_dict)#, style='size_code')
    plt.legend(title='Size code')
    ax.legend([],[], frameon=False)
    #sns.despine()

    plt.ylabel('Relative abundance', fontsize=12)
    plt.xlabel('Time (week)', fontsize=12)


    sns.set_style("ticks")
    plt.subplot(224)
    ax = sns.lineplot(data = sfd_f, x='weekn', y='ranktot', hue='size_code', palette=palette_dict)#, style='size_code')
    plt.legend(title='Size code')
    ax.legend([],[], frameon=False)
    plt.ylim(reversed(plt.ylim()))
    #sns.despine()

    plt.ylabel('Rank', fontsize=12)
    plt.xlabel('Time (week)', fontsize=12)
    plt.tight_layout()

    plt.savefig('outputs/'+comm+'/D'+str(depth)+'_lineplot'+ fid +'allplots.png', dpi=200, bbox_inches="tight")


def save_individual_plots(comm,depth,fid,newseparated):

    sizecodes = ['S', 'L', 'W', 'SL']
    palette_colors = sns.color_palette()
    palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}

    sfd=newseparated[newseparated.depth==depth]
    sfd['weekfid'] = sfd["weekn"].astype(str) + sfd["feature_id"].astype(str)
    sfd['avg_p_id'] = sfd['ratio'].groupby(sfd['weekfid']).transform('mean')
    sfd['diff_p_id'] = sfd['ratio'] - sfd['avg_p_id']

    sfd_f=sfd[sfd.feature_id==fid]

    tax = sfd.loc[sfd["feature_id"] == fid].iloc[0]["Taxon"]
    rk = sfd.loc[sfd["feature_id"] == fid].iloc[0]["rank"]

    ttl = sfd_f['Taxon'].iloc[0]

    sns.set(rc={"figure.figsize":(7, 4)})
    sns.set_style("ticks")
    ax=sns.lineplot(data=sfd_f, x="weekn", y="diff_p_id", hue="size_code", palette=palette_dict)#, hue="size_code")
    plt.legend(title='Size code')
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), frameon=False)

    plt.title(ttl)
    plt.ylabel('ΔRA per week')
    plt.xlabel('Week number')
    plt.savefig('outputs/'+comm+'/D'+str(depth)+fid+'avg_per_week.png', dpi=200, bbox_inches="tight")
    plt.clf()
    plt.cla()
    plt.close()


    sfd_f['SCfid'] = sfd_f["size_code"].astype(str) + sfd_f["feature_id"].astype(str)
    sfd_f['avg_p_sc'] = sfd_f['ratio'].groupby(sfd_f['SCfid']).transform('mean')
    sfd_f['diff_p_sc'] = sfd_f['ratio'] - sfd_f['avg_p_sc']


    sns.set(rc={"figure.figsize":(7, 4)})
    sns.set_style("ticks")
    ax=sns.lineplot(data=sfd_f, x="weekn", y="diff_p_sc", hue="size_code", palette=palette_dict)#, hue="size_code")
    plt.legend(title='Size code')
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), frameon=False)

    plt.title(ttl)
    plt.ylabel('ΔRA per size fraction')
    plt.xlabel('Week number')
    plt.savefig('outputs/'+comm+'/D'+str(depth)+fid+'_avg_per_SC.png', dpi=200, bbox_inches="tight")
    plt.clf()
    plt.cla()
    plt.close()


    sns.set(rc={"figure.figsize":(7, 4)})
    sns.set_style("ticks")
    ax = sns.lineplot(data = sfd_f, x='weekn', y='ratio', hue='size_code', palette=palette_dict)#, style='size_code')
    plt.legend(title='Size code')
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), frameon=False)
    #sns.despine()

    plt.ylabel('Relative abundance', fontsize=12)
    plt.xlabel('Time (week)', fontsize=12)
    plt.title('Relative abundance of '+ tax + '('+ str(rk) +')')

    plt.savefig('outputs/'+comm+'/D'+str(depth)+'_lineplot'+ tax +'RA.png', dpi=200, bbox_inches="tight")
    plt.clf()
    plt.cla()
    plt.close()


    sns.set(rc={"figure.figsize":(7, 4)})
    sns.set_style("ticks")
    ax = sns.lineplot(data = sfd_f, x='weekn', y='ranktot', hue='size_code', palette=palette_dict)#, style='size_code')
    plt.legend(title='Size code')
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), frameon=False)
    plt.ylim(reversed(plt.ylim()))
    #sns.despine()

    plt.ylabel('Rank', fontsize=12)
    plt.xlabel('Time (week)', fontsize=12)
    plt.title('Relative abundance of '+ tax + '('+ str(rk) +')')

    plt.savefig('outputs/'+comm+'/D'+str(depth)+'_lineplot'+ tax +'RANK.png', dpi=200, bbox_inches="tight")


def boxplot_depth(separated, comm, depth, ycolumn, yaxislabel='def'):

    if yaxislabel != 'def':
        ycol = ycolumn

    if depth == 'all':
        sfd = separated.copy()
    else:
        sfd=separated[separated.depth==depth]

    #sfd_S = sfd[['size_code', 'nASVs', 'weekn']].copy()
    #sfd_S = sfd_S.drop_duplicates()
    #sdfpv = sfd_S.pivot(index='weekn', columns='size_code', values='nASVs')
    #fvalue, pvalue = stats.f_oneway(sdfpv['L'], sdfpv['S'], sdfpv['W'])

    sfd_LM = sfd[['size_code', 'nASVs']].copy()
    sfd_LM = sfd_LM.drop_duplicates()
    print(sfd_LM['size_code'].unique())
    lm = sfa.ols('nASVs ~ C(size_code)', data=sfd_LM).fit()
    anova = sa.stats.anova_lm(lm)
    results = spPH.posthoc_ttest(sfd_LM, val_col='nASVs', group_col='size_code', p_adjust='holm')

    if 'SL' in separated['size_code'].unique():
        #sizecode palette codes
        sizecodes = ['S', 'L', 'W', 'SL', 'P']
        palette_colors = sns.color_palette()
        palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}

    else:
        #define color palettes
        sizecodes = ['S', 'L', 'W', 'P']
        palette_colors = sns.color_palette()
        palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}



    #plot
    sns.set(rc={"figure.figsize":(4, 3)})
    sns.set_style("ticks")
    sns.boxplot(data=sfd, x="size_code", y=ycolumn, palette=palette_dict, order=sizecodes)#, hue="size_code")
    sns.despine()
    plt.ylabel(yaxislabel, fontsize=20)
    plt.xlabel('Size fraction', fontsize=20)

    #g.tick_params(labelsize=15)
    plt.savefig('outputs/'+comm+'/D'+str(depth)+'_adboxplot.png', dpi=200, bbox_inches="tight")
    plt.clf()
    plt.cla()
    plt.close()


    sns.set(rc={"figure.figsize":(7, 3)})
    sns.set_style("ticks")
    ax=sns.barplot(data=sfd, x="weekn", y="diff", hue="size_code", palette=palette_dict,
                  capsize=.15, errwidth=0.5)#, hue="size_code")
    sns.despine()
    plt.ylabel('Number of ASVs relative to weekly average', fontsize=20)
    plt.xlabel('Week number', fontsize=20)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.savefig('outputs/'+comm+'/D'+str(depth)+'_avgbarplot.png', dpi=200, bbox_inches="tight")

    plt.clf()


#     glue = sfd[['size_code', 'weekn', 'diff']].copy()
#     glue = glue.drop_duplicates()
#     glue = glue.pivot(index="size_code", columns="weekn", values="diff")
#     floored_data = glue.apply(np.floor)
#     sns.set_style('ticks')
#     plt.figure(figsize=(8, 2))
#     cmap = sns.diverging_palette(240,240, as_cmap=True)
#     ax = sns.heatmap(floored_data, yticklabels=True, linewidths=.5, annot=True, annot_kws={"fontsize":8},
#                     cmap = cmap)
#     plt.savefig('outputs/'+comm_id+'/heatmap_nasv_change_d'+str(depth)+'_annot.png', bbox_inches='tight', dpi=300)
#     plt.clf()

#     ax = sns.heatmap(floored_data, fmt='.1f', yticklabels=True, linewidths=.5,
#                     cmap = cmap)
#     plt.savefig('outputs/'+comm_id+'/heatmap_nasv_change_d'+str(depth)+'.png', bbox_inches='tight', dpi=300)


    plt.clf()
    sns.set(rc={"figure.figsize":(7, 3)})
    sns.set_style("ticks")
    ax=sns.lineplot(x = "weekn", y = ycolumn, data=sfd, hue="size_code", palette=palette_dict)
    sns.despine()
    plt.ylabel(yaxislabel, fontsize=20)
    plt.xlabel('Week', fontsize=20)
    plt.legend(title='Size fraction')
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.savefig('outputs/'+comm+'/D'+str(depth)+'_adlineplot.png', dpi=200, bbox_inches="tight")

    return anova, results


# In[ ]:

def subplots(comm, depth, fids, pltyp, newseparated):
    sizecodes = ['S', 'L', 'W', 'SL']
    palette_colors = sns.color_palette()
    palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}

    fig, ax = plt.subplots(len(fids), 1, sharex=True)
    plt.subplots_adjust(wspace=0, hspace=0)

    d1 = newseparated[newseparated.depth == depth]

    sns.set(rc={"figure.figsize":(4, 7)})
    sns.set_style("ticks")


    for i in range(len(fids)):
        d1_fid = d1[d1.feature_id == fids[i]]
        if comm == 'chloroplast':
            ttl = d1.loc[d1["feature_id"] == fids[i]].iloc[0]["PRSpecies"]
        else:
            ttl = d1.loc[d1["feature_id"] == fids[i]].iloc[0]["Genus"]

        d1_fid['SCfid'] = d1_fid["size_code"].astype(str) + d1_fid["feature_id"].astype(str)
        d1_fid['avg_p_sc'] = d1_fid['ratio'].groupby(d1_fid['SCfid']).transform('mean')
        d1_fid['diff_p_sc'] = d1_fid['ratio'] - d1_fid['avg_p_sc']

        d1_fid['weekfid'] = d1_fid["weekn"].astype(str) + d1_fid["feature_id"].astype(str)
        d1_fid['avg_p_id'] = d1_fid['ratio'].groupby(d1_fid['weekfid']).transform('mean')
        d1_fid['diff_p_id'] = d1_fid['ratio'] - d1_fid['avg_p_id']

        #pick plot type
        if pltyp == 1:
            ax1=sns.lineplot(data=d1_fid, x="weekn", y="diff_p_sc", hue="size_code", palette=palette_dict, ax=ax[i])
        if pltyp == 2:
            ax1=sns.lineplot(data=d1_fid, x="weekn", y="diff_p_id", hue="size_code", palette=palette_dict, ax=ax[i])
        if pltyp == 3:
            ax1 = sns.lineplot(data = d1_fid, x='weekn', y='ratio', hue='size_code', palette=palette_dict, ax=ax[i])
        if pltyp == 4:
            ax1 = sns.lineplot(data = d1_fid, x='weekn', y='ranktot', hue='size_code', palette=palette_dict, ax=ax[i]) #

        ax1.tick_params(bottom=False)
        ax1.set(ylabel=None)
        ax1.get_legend().remove()
        ax1.set_ylabel(ttl, rotation=0, labelpad=85)

        for _,s in ax1.spines.items():
            s.set_linewidth(0.5)
            s.set_color('grey')

    plt.savefig('outputs/'+comm+'/D'+str(depth)+'_lineplots'+ str(pltyp) +'.png', dpi=200, bbox_inches="tight")


def upsetprep(comm, level, separated, depth):

    depths = [1, 5, 10, 30, 60]

    cumulab = separated[['feature_frequency', 'depth', 'size_code', level]].copy()
    cumulab1 = cumulab.groupby([level]).agg({'feature_frequency':sum})

    if depth != 'all':
        sfd=separated[separated.depth==depth]
    else:
        sfd=separated.copy()

    toptaxa = sfd[['feature_frequency', 'Taxon', 'size_code', 'depth','weekn', level]].copy()
    toptaxa = toptaxa.drop_duplicates()
    df_agg = toptaxa.groupby(['size_code',level, 'depth']).agg({'feature_frequency':sum})

    df_agg = df_agg.reset_index()
    df_agg['set_name'] = df_agg['size_code']+df_agg['depth'].astype(str)

    resultpivot = df_agg.pivot_table(index=level, columns='set_name', values='feature_frequency')
    resultpivot = resultpivot.fillna(0)
    resultpivot[resultpivot != 0] = 1
    tosave = pd.merge(resultpivot, cumulab1, left_index=True, right_index=True)
    tosave.to_csv('csvs/'+comm+'/'+level+'_d'+str(depth)+'_relab.csv')


    #make json
    data = {
        "file": "https://raw.githubusercontent.com/dianahaider/size_fractions/main/csvs/"+comm+'/'+level+'_d'+str(depth)+'_relab.csv',
        "name": comm + level,
        "header": 0,
        "separator": ",",
        "skip": 0,
        "meta":[
            {"type":"id", "index":0, "name":"Name"},
            {"type":"integer", "index":4, "name":"Rel. ab."}
        ],
        "sets": [
            {"format": "binary", "start":1, "end": 3}
        ]
    }

    with open('json/'+comm+'/'+level+'_d'+str(depth)+'.json', 'w') as f:
        json.dump(data, f)


# In[ ]:


def plot_per_fid(comm, separated, depth, fid):


    if 'SL' in separated['size_code'].unique():
        #sizecode palette codes
        sizecodes = ['S', 'L', 'W', 'SL']
        palette_colors = sns.color_palette()
        palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}

    else:
        #sizecode palette codes
        sizecodes = ['S', 'L', 'W']
        palette_colors = sns.color_palette()
        palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}

    sfd=separated[separated.depth==depth]
    sfd['weekfid'] = sfd["weekn"].astype(str) + sfd["feature_id"].astype(str)
    sfd['avg_p_id'] = sfd['ratio'].groupby(sfd['weekfid']).transform('mean')
    sfd['diff_p_id'] = sfd['ratio'] - sfd['avg_p_id']

    sfd_f=sfd[sfd.feature_id==fid]

    ttl = sfd_f['Taxon'].iloc[0]

    sns.set(rc={"figure.figsize":(7, 3)})
    ax=sns.lineplot(data=sfd_f, x="weekn", y="diff_p_id", hue="size_code", palette=palette_dict)#, hue="size_code")
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.title(ttl)
    plt.ylabel('Ratio difference')
    plt.xlabel('Week number')
    plt.savefig('outputs/'+comm+'/D'+str(depth)+fid+'.png', dpi=200, bbox_inches="tight")


def run_ancom(comm, separated, sfdclr, depth, ancomcol):
    # Assign group labels to the samples based on depth
    sfd = separated[separated.depth == depth]
    df_ancom = sfd[['sampleid', ancomcol]].copy().drop_duplicates().set_index('sampleid')

    # Run ANCOM on the CLR-transformed table
    results = ancom(table=sfdclr, grouping=df_ancom[ancomcol], multiple_comparisons_correction='holm-bonferroni')
    DAresults = results[0].copy()  # Contains W values and significance

    # Taxonomy information for chloroplasts or other communities
    if comm == 'chloroplast':
        taxonomy = sfd[['feature_id', 'PRConfidence', 'PRTaxon', 'PRSpecies']].drop_duplicates()
    else:
        taxonomy = sfd[['feature_id', 'Confidence', 'Taxon', 'Phylum', 'Class', 'Family', 'Genus', 'Species']].drop_duplicates()
    DAtaxonomy = pd.merge(DAresults, taxonomy, on="feature_id", how="left")

    # Extract the significant features for further analysis
    DARejected_SC = DAtaxonomy.loc[DAtaxonomy['Reject null hypothesis'] == True].sort_values(by=['W'])

    # Calculate CLR mean values for each feature within each group
    clr_means = sfdclr.copy()
    clr_means['group'] = df_ancom[ancomcol]  # Assign group labels
    clr_means_grouped = clr_means.groupby('group').mean().T  # Mean CLR per group

    # Calculate the standard deviation of CLR mean values across groups
    clr_means_grouped['clr_mean_std_dev'] = clr_means_grouped.std(axis=1)

    # Merge the W values and CLR mean standard deviations for the volcano plot
    volcano_df = DAresults.copy()
    volcano_df['clr_mean_std_dev'] = clr_means_grouped['clr_mean_std_dev']
    volcano_df = pd.merge(volcano_df, taxonomy, on="feature_id", how="left")  # Add taxonomy for labeling

    # Define threshold for highlighting significant features
    w_threshold = volcano_df['W'].quantile(0.95)  # e.g., 95th percentile of W as significance threshold
    std_dev_threshold = 0.5  # Example threshold for CLR standard deviation



    # Create the Volcano Plot with Enhanced Visualization
    plt.figure(figsize=(12, 8))
    sns.scatterplot(
        x='clr_mean_std_dev',
        y='W',
        data=volcano_df,
        hue='Reject null hypothesis',
        palette={True: "red", False: "blue"},
        alpha=0.6,
        edgecolor="k",
        s=80
    )

    # Highlight top significant taxa
    for _, row in volcano_df[(volcano_df['W'] > w_threshold) |
                             (volcano_df['clr_mean_std_dev'] > std_dev_threshold)].iterrows():
        # Determine the label based on the community type
        if comm == 'chloroplast' and pd.notnull(row['PRSpecies']):
            label = row['PRSpecies']
        elif pd.notnull(row['Genus']):
            label = row['Genus']
        else:
            label = row['feature_id']

        plt.text(
            row['clr_mean_std_dev'],
            row['W'],
            label,
            fontsize=9,
            ha='right' if row['clr_mean_std_dev'] > 0 else 'left',
            color="black"
        )

    # Add plot labels and title
    plt.axhline(y=w_threshold, color='grey', linestyle='--', label='Significance Threshold (95th percentile)')
    plt.xlabel('Standard Deviation of CLR Mean Abundance')
    plt.ylabel('W Statistic')
    plt.title('Enhanced ANCOM Volcano Plot for Multiple Groups')
    plt.legend(title="Significance")
    plt.grid(True)
    plt.tight_layout()

    # Save the plot as an image file
    plt.savefig('outputs/'+comm+'/D'+str(depth)+'volcano_plot.png', dpi=300, bbox_inches="tight")
    plt.close()  # Close the plot to free memory if running in a loop or batch processing

    prcentile = results[1].copy()  # Contains percentile abundance per group

    return DAtaxonomy, DARejected_SC, prcentile


def run_ancoms(comm, separated, sfdclr, depth, ancomcol, threshold=0):
    #create a new df where we remove the low abundance ASVs optionally
    news2 = newseparated.drop(newseparated[(newseparated['ratio'] <threshold )].index)

    #re-calculate ratios when removing chloroplast for 16S
    news2['Total'] = news2['feature_frequency'].groupby(news2['sampleid']).transform('sum')
    news2['ratio'] = news2['feature_frequency']/news2['Total']

    #sfd is used to assign group labels to the samples
    sfd=separated[separated.depth==depth]
    df_ancom = sfd[['sampleid', ancomcol]].copy()
    df_ancom = df_ancom.drop_duplicates()
    df_ancom = df_ancom.set_index('sampleid')

    #sfdclr is the clr transformed table (columns are asvs, rows are samples)
    results = ancom(table=sfdclr, grouping=df_ancom[ancomcol], multiple_comparisons_correction='holm-bonferroni')
    DAresults = results[0].copy()

    if comm == 'chloroplast':
        taxonomy = sfd[['feature_id', 'PRConfidence', 'PRTaxon', 'PRSpecies']].copy()
        taxonomy = taxonomy.drop_duplicates()
    else:
        taxonomy = sfd[['feature_id', 'Confidence', 'Taxon', 'Phylum', 'Class', 'Family', 'Genus', 'Species']].copy()
        taxonomy = taxonomy.drop_duplicates()
    DAtaxonomy = pd.merge(DAresults, taxonomy, on="feature_id", how="left")

    DARejected_SC = DAtaxonomy.loc[DAtaxonomy['Reject null hypothesis'] == True]
    DARejected_SC.sort_values(by=['W'])


    prcentile = results[1].copy()

    return DAtaxonomy, DARejected_SC, prcentile


# In[ ]:


subtitle = 'From Jan7 2022 to Apr27 2022'
def plot_stackedbar_(comm, df, labels, colors, title, subtitle, level):

    fields = df.columns.tolist()

    # figure and axis
    fig, ax = plt.subplots(1, figsize=(8, 5))
# plot bars
    left = len(df) * [0]
    for idx, name in enumerate(fields):
        plt.barh(df.index, df[name], left = left, color=colors[idx])
        left = left + df[name]
# title and subtitle
    plt.title(title, loc='left')
    #plt.text(0, ax.get_yticks()[-1] + 0.75, subtitle)
# legend
    plt.legend(labels, bbox_to_anchor=([1, 1, 0, 0]), ncol=1, frameon=False)
# remove spines
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
# format x ticks
    xticks = np.arange(0,110,10)
    xlabels = ['{}%'.format(i) for i in np.arange(0,101,10)]
    plt.xticks(xticks, xlabels)
# adjust limits and draw grid lines
    plt.ylim(-0.5, ax.get_yticks()[-1] + 0.5)
    ax.xaxis.grid(color='gray', linestyle='dashed')
    plt.gca().invert_yaxis()
    plt.ylabel("Depth (m)")
    plt.savefig('outputs/'+comm+'/'+level+'alldepths_stacked_perc_weighted.png', dpi=200, bbox_inches="tight")
    plt.show()


# In[ ]:


subtitle = 'From Jan7 2022 to Apr27 2022'
def plot_stackedbar_p(comm, df, labels, colors, title, subtitle, level):

    fields = df.columns.tolist()

    # figure and axis
    fig, ax = plt.subplots(1, figsize=(8, 5))
# plot bars
    left = len(df) * [0]
    for idx, name in enumerate(fields):
        plt.barh(df.index, df[name], left = left, color=colors[idx])
        left = left + df[name]
# title and subtitle
    plt.title(title, loc='left')
    #plt.text(0, ax.get_yticks()[-1] + 0.75, subtitle)
# legend
    plt.legend(labels, bbox_to_anchor=([1, 1, 0, 0]), ncol=1, frameon=False)
# remove spines
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
# format x ticks
    xticks = np.arange(0,110,10)
    xlabels = ['{}%'.format(i) for i in np.arange(0,101,10)]
    plt.xticks(xticks, xlabels)
# adjust limits and draw grid lines
    plt.ylim(-0.5, ax.get_yticks()[-1] + 0.5)
    ax.xaxis.grid(color='gray', linestyle='dashed')
    plt.gca().invert_yaxis()
    plt.ylabel("Depth (m)")
    plt.savefig('outputs/'+comm+'/'+level+'alldepths_stacked_perc_weighted.png', dpi=200, bbox_inches="tight")
    plt.show()


# In[ ]:


subtitle = 'From Jan7 2022 to Apr27 2022'
def plot_stackedbar_p(comm, df, labels, colors, title, subtitle, level):

    fields = df.columns.tolist()

    # figure and axis
    fig, ax = plt.subplots(1, figsize=(8, 5))
# plot bars
    left = len(df) * [0]
    for idx, name in enumerate(fields):
        plt.barh(df.index, df[name], left = left, color=colors[idx])
        left = left + df[name]
# title and subtitle
    plt.title(title, loc='left')
    #plt.text(0, ax.get_yticks()[-1] + 0.75, subtitle)
# legend
    plt.legend(labels, bbox_to_anchor=([1, 1, 0, 0]), ncol=1, frameon=False)
# remove spines
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
# format x ticks
    xticks = np.arange(0,110,10)
    xlabels = ['{}%'.format(i) for i in np.arange(0,101,10)]
    plt.xticks(xticks, xlabels)
# adjust limits and draw grid lines
    plt.ylim(-0.5, ax.get_yticks()[-1] + 0.5)
    ax.xaxis.grid(color='gray', linestyle='dashed')
    plt.gca().invert_yaxis()
    plt.ylabel("Depth (m)")
    plt.savefig('outputs/'+comm+'/'+level+'alldepths_stacked_perc_weighted.png', dpi=200, bbox_inches="tight")
    plt.show()


# In[ ]:


def calcperc_defrac(comm, separated, level, dfplot_unweighted):


    depths = [1, 5, 10, 30, 60]

    level = level

    #make an empty df to plot
    dfplot = pd.DataFrame(columns=['Depth', 'SF', 'NSF', 'DFr', 'Both'])

    for d in range(len(depths)):
        sfd=separated[separated.depth==depths[d]]
        toptaxa = sfd[[level, 'feature_frequency', 'Taxon', 'size_code', 'weekn']].copy()

        toptaxa = toptaxa.drop_duplicates()
        df_agg = toptaxa.groupby(['size_code',level]).agg({'feature_frequency':sum})

        df_agg = df_agg.reset_index()
        resultpivot = df_agg.pivot_table(index=level, columns='size_code', values='feature_frequency')
        resultpivot = resultpivot.fillna(0)

        df1 = resultpivot.copy()

        df = resultpivot[['L', 'S', 'W']].copy()
        Sonly = df[(df['L'] == 0) & (df['W'] == 0)]
        Wonly = df[(df['L'] == 0) & (df['S'] == 0)]
        Lonly = df[(df['S'] == 0) & (df['W'] == 0)]
        LW = df[(df['S'] == 0) & (df['W'] != 0) & (df['L'] != 0)]
        LS = df[(df['W'] == 0) & (df['S'] != 0) & (df['L'] != 0)]
        SW = df[(df['W'] != 0) & (df['S'] != 0) & (df['L'] == 0)]
        LSW = df[~(df == 0).any(axis=1)]

        DFr = df1[(df1['SL'] != 0)]
        DFr = DFr[['SL']].copy()

        total = resultpivot.to_numpy().sum()

        SFdf = Lonly, LS, Sonly
        SF = pd.concat(SFdf)
        SF_value = SF.to_numpy().sum()/total *100

        Sonly_value = Sonly.to_numpy().sum()/total *100


        Bothdf = LW, LSW, SW
        Both = pd.concat(Bothdf)
        Both_value = Both.to_numpy().sum()/total *100

        Wonly_value = Wonly.to_numpy().sum()/total *100

        Sonly_value = Sonly.to_numpy().sum()/total *100
        Lonly_value = Lonly.to_numpy().sum()/total *100
        LS_value = LS.to_numpy().sum()/total *100
        DFr_value = DFr.to_numpy().sum()/total *100

        dfplot.loc[d,'Depth'] = depths[d]
        dfplot.loc[d,'SF'] = SF_value
        dfplot.loc[d,'NSF'] = Wonly_value
        dfplot.loc[d,'DFr'] = DFr_value
        dfplot.loc[d,'Both'] = Both_value

        dfplot_unweighted.loc[d,'Depth'] = depths[d]
        dfplot_unweighted.loc[d,'SF'] = len(Lonly) + len(Sonly) + len(LS)
        dfplot_unweighted.loc[d,'NSF'] = len(Wonly)
        dfplot_unweighted.loc[d,'Both'] = len(LW) + len(SW) + len(LSW)


        venn3(subsets = (len(Lonly), len(Sonly), len(LS), len(Wonly), len(LW), len(SW), len(LSW)), set_labels = ('Large >3μm', 'Small 3-02μm', 'Whole water <0.22μm'), alpha = 0.5);

        plt.savefig("outputs/"+comm+"/D"+str(depths[d])+level+"_venn.png")
        plt.clf()
        plt.cla()
        plt.close()

    dfplot['Depth'] = dfplot['Depth'].astype(str)
    dfplot = dfplot.set_index('Depth')

    return dfplot, level, dfplot_unweighted


# In[ ]:


def calcperc_defrac_unweighted(comm, separated, level):
    depths = [1, 5, 10, 30, 60]

    level = level

    dfplot = pd.DataFrame(columns=['Depth', 'SF', 'NSF', 'Both'])
    dfplot_unweighted = pd.DataFrame(columns=['Depth', 'SF', 'NSF', 'Both'])

    for d in range(len(depths)):
        sfd=separated[separated.depth==depths[d]]
        toptaxa = sfd[[level, 'feature_frequency', 'Taxon', 'size_code', 'weekn']].copy()

        toptaxa = toptaxa.drop_duplicates()
        df_agg = toptaxa.groupby(['size_code',level]).agg({'feature_frequency':sum})

        df_agg = df_agg.reset_index()
        resultpivot = df_agg.pivot_table(index=level, columns='size_code', values='feature_frequency')
        resultpivot = resultpivot.fillna(0)

        df1 = resultpivot.copy()

        df = resultpivot[['L', 'S', 'W']].copy()
        Sonly = df[(df['L'] == 0) & (df['W'] == 0)]
        Wonly = df[(df['L'] == 0) & (df['S'] == 0)]
        Lonly = df[(df['S'] == 0) & (df['W'] == 0)]
        LW = df[(df['S'] == 0) & (df['W'] != 0) & (df['L'] != 0)]
        LS = df[(df['W'] == 0) & (df['S'] != 0) & (df['L'] != 0)]
        SW = df[(df['W'] != 0) & (df['S'] != 0) & (df['L'] == 0)]
        LSW = df[~(df == 0).any(axis=1)]

        total = len(resultpivot)

        SFdf = Lonly, LS, Sonly
        SF = pd.concat(SFdf)
        SF_value = len(SF)/total *100

        Sonly_value = len(Sonly)/total *100


        Bothdf = LW, LSW, SW
        Both = pd.concat(Bothdf)
        Both_value = len(Both)/total *100

        Wonly_value = len(Wonly)/total *100

        Sonly_value = len(Sonly)/total *100
        Lonly_value = len(Lonly)/total *100
        LS_value = len(LS)/total *100

        dfplot.loc[d,'Depth'] = depths[d]
        dfplot.loc[d,'SF'] = SF_value
        dfplot.loc[d,'NSF'] = Wonly_value
        dfplot.loc[d,'Both'] = Both_value

        dfplot_unweighted.loc[d,'Depth'] = depths[d]
        dfplot_unweighted.loc[d,'SF'] = Lonly_value + Sonly_value + LS_value
        dfplot_unweighted.loc[d,'NSF'] = Wonly_value
        dfplot_unweighted.loc[d,'Both'] = 100 - (Lonly_value + Sonly_value + LS_value + Wonly_value)


        venn3(subsets = (len(Lonly), len(Sonly), len(LS), len(Wonly), len(LW), len(SW), len(LSW)), set_labels = ('Large >3μm', 'Small 3-02μm', 'Whole water <0.22μm'), alpha = 0.5);

        plt.savefig("outputs/"+comm+"/D"+str(depths[d])+level+"_venn.png")
        plt.clf()
        plt.cla()
        plt.close()

    dfplot['Depth'] = dfplot['Depth'].astype(str)
    dfplot = dfplot.set_index('Depth')

    return dfplot, level, dfplot_unweighted


# In[ ]:


def calcperc(comm, separated, level):

    depths = [1, 5, 10, 30, 60]

    level = level

    dfplot = pd.DataFrame(columns=['Depth', 'Sonly', 'Lonly', 'LS', 'NSF'])

    #save venn of all depths together
    toptaxa = separated[[level, 'feature_frequency', 'Taxon', 'size_code', 'weekn']].copy()

    toptaxa = toptaxa.drop_duplicates()
    df_agg = toptaxa.groupby(['size_code',level]).agg({'feature_frequency':sum})

    df_agg = df_agg.reset_index()
    resultpivot = df_agg.pivot_table(index=level, columns='size_code', values='feature_frequency')
    resultpivot = resultpivot.fillna(0)

    df = resultpivot.copy()

    Sonly = df[(df['L'] == 0) & (df['W'] == 0)]
    Wonly = df[(df['L'] == 0) & (df['S'] == 0)]
    Lonly = df[(df['S'] == 0) & (df['W'] == 0)]
    LW = df[(df['S'] == 0) & (df['W'] != 0) & (df['L'] != 0)]
    LS = df[(df['W'] == 0) & (df['S'] != 0) & (df['L'] != 0)]
    SW = df[(df['W'] != 0) & (df['S'] != 0) & (df['L'] == 0)]
    LSW = df[~(df == 0).any(axis=1)]

    D_taxlist = pd.DataFrame.from_dict({'Sonly': Sonly.index.tolist(), 'Wonly': Wonly.index.tolist(), 'Lonly': Lonly.index.tolist(),
            'LW': LW.index.tolist(), 'LS':LS.index.tolist(), 'SW':SW.index.tolist(), 'LSW':LSW.index.tolist()}, orient='index').T

    D_taxlist.to_csv('outputs/'+comm+'/all_taxlist.csv')


    total = df.to_numpy().sum()

    Wonly_value = Wonly.to_numpy().sum()/total *100

    Sonly_value = Sonly.to_numpy().sum()/total *100
    Lonly_value = Lonly.to_numpy().sum()/total *100
    LS_value = LS.to_numpy().sum()/total *100

    LW_value = LW.to_numpy().sum()/total *100
    SW_value = SW.to_numpy().sum()/total *100
    LSW_value = LSW.to_numpy().sum()/total *100

    vd = venn3(subsets = (len(Lonly), len(Sonly), len(LS), len(Wonly), len(LW), len(SW), len(LSW)),
        set_labels = ('Large >3μm', 'Small 0.2-3μm', 'Whole water >0.2μm'),
        set_colors=("#d55e00", "#5975a4", "#5f9e6e"),
        alpha=1);
    for text in vd.set_labels:
        text.set_fontsize(16)
    for text in vd.subset_labels:
        text.set_fontsize(18)

    vd2 = venn3(subsets = (Lonly_value, Sonly_value, LS_value, Wonly_value, LW_value, SW_value, LSW_value),
        set_labels = ('Large >3μm', 'Small 3-02μm', 'Whole water <0.22μm'),
        set_colors=("#d55e00", "#5975a4", "#5f9e6e"),
        alpha=1);

    plt.savefig("outputs/"+comm+"/alldepths_"+level+"_venn_2.png")
    plt.clf()
    plt.cla()
    plt.close()


    for d in range(len(depths)):


        sfd=separated[separated.depth==depths[d]]
        toptaxa = sfd[[level, 'feature_frequency', 'Taxon', 'size_code', 'weekn']].copy()

        toptaxa = toptaxa.drop_duplicates()
        df_agg = toptaxa.groupby(['size_code',level]).agg({'feature_frequency':sum})

        df_agg = df_agg.reset_index()
        resultpivot = df_agg.pivot_table(index=level, columns='size_code', values='feature_frequency')
        resultpivot = resultpivot.fillna(0)

        df = resultpivot.copy()

        Sonly = df[(df['L'] == 0) & (df['W'] == 0)]
        Wonly = df[(df['L'] == 0) & (df['S'] == 0)]
        Lonly = df[(df['S'] == 0) & (df['W'] == 0)]
        LW = df[(df['S'] == 0) & (df['W'] != 0) & (df['L'] != 0)]
        LS = df[(df['W'] == 0) & (df['S'] != 0) & (df['L'] != 0)]
        SW = df[(df['W'] != 0) & (df['S'] != 0) & (df['L'] == 0)]
        LSW = df[~(df == 0).any(axis=1)]

        D_taxlist = pd.DataFrame.from_dict({'Sonly': Sonly.index.tolist(), 'Wonly': Wonly.index.tolist(), 'Lonly': Lonly.index.tolist(),
            'LW': LW.index.tolist(), 'LS':LS.index.tolist(), 'SW':SW.index.tolist(), 'LSW':LSW.index.tolist()}, orient='index').T

        D_taxlist.to_csv('outputs/'+comm+'/'+str(depths[d])+'_taxlist.csv')


        total = df.to_numpy().sum()

        SFdf = Lonly, LS, Sonly
        SF = pd.concat(SFdf)
        SF_value = SF.to_numpy().sum()/total *100

        Sonly_value = Sonly.to_numpy().sum()/total *100


        Bothdf = LW, LSW, SW
        Both = pd.concat(Bothdf)
        Both_value = Both.to_numpy().sum()/total *100

        Wonly_value = Wonly.to_numpy().sum()/total *100

        Sonly_value = Sonly.to_numpy().sum()/total *100
        Lonly_value = Lonly.to_numpy().sum()/total *100
        LS_value = LS.to_numpy().sum()/total *100

        LW_value = LW.to_numpy().sum()/total *100
        SW_value = SW.to_numpy().sum()/total *100
        LSW_value = LSW.to_numpy().sum()/total *100


        dfplot.loc[d,'Depth'] = depths[d]
        dfplot.loc[d,'Sonly'] = Sonly_value
        dfplot.loc[d,'Lonly'] = Lonly_value
        dfplot.loc[d,'LS'] = LS_value
        dfplot.loc[d,'NSF'] = Wonly_value
        dfplot.loc[d,'Both'] = Both_value

        vd3 = venn3(subsets = (Lonly_value, Sonly_value, LS_value, Wonly_value, LW_value, SW_value, LSW_value),
            set_labels = ('Large >3μm', 'Small 3-02μm', 'Whole water <0.22μm'),
            set_colors=("#d55e00", "#5975a4", "#5f9e6e"),
            alpha=1);


        plt.savefig("outputs/"+comm+"/D"+str(depths[d])+level+"_venn_2.png")
        plt.clf()
        plt.cla()
        plt.close()


    dfplot['Depth'] = dfplot['Depth'].astype(str)
    dfplot = dfplot.set_index('Depth')

    return dfplot, level


# In[ ]:


def calcperc_SLNSF(comm, separated, level):

    depths = [1, 5, 10, 30, 60]

    level = level

    dfplot = pd.DataFrame(columns=['Depth', 'Sonly', 'Lonly', 'LS', 'NSF'])

    for d in range(len(depths)):
        sfd=separated[separated.depth==depths[d]]
        toptaxa = sfd[[level, 'feature_frequency', 'Taxon', 'size_code', 'weekn']].copy()

        toptaxa = toptaxa.drop_duplicates()
        df_agg = toptaxa.groupby(['size_code',level]).agg({'feature_frequency':sum})

        df_agg = df_agg.reset_index()
        resultpivot = df_agg.pivot_table(index=level, columns='size_code', values='feature_frequency')
        resultpivot = resultpivot.fillna(0)

        df = resultpivot.copy()

        Sonly = df[(df['L'] == 0) & (df['W'] == 0)]
        Wonly = df[(df['L'] == 0) & (df['S'] == 0)]
        Lonly = df[(df['S'] == 0) & (df['W'] == 0)]
        LW = df[(df['S'] == 0) & (df['W'] != 0) & (df['L'] != 0)]
        LS = df[(df['W'] == 0) & (df['S'] != 0) & (df['L'] != 0)]
        SW = df[(df['W'] != 0) & (df['S'] != 0) & (df['L'] == 0)]
        LSW = df[~(df == 0).any(axis=1)]

        total = df.to_numpy().sum()

        SFdf = Lonly, LS, Sonly
        SF = pd.concat(SFdf)
        SF_value = SF.to_numpy().sum()/total *100

        Sonly_value = Sonly.to_numpy().sum()/total *100


        Bothdf = LW, LSW, SW
        Both = pd.concat(Bothdf)
        Both_value = Both.to_numpy().sum()/total *100

        Wonly_value = Wonly.to_numpy().sum()/total *100

        Sonly_value = Sonly.to_numpy().sum()/total *100
        Lonly_value = Lonly.to_numpy().sum()/total *100
        LS_value = LS.to_numpy().sum()/total *100

        NewTotal = Sonly_value + Lonly_value + LS_value + Wonly_value


        dfplot.loc[d,'Depth'] = depths[d]
        dfplot.loc[d,'Sonly'] = Sonly_value
        dfplot.loc[d,'Lonly'] = Lonly_value
        dfplot.loc[d,'LS'] = LS_value
        dfplot.loc[d,'NSF'] = Wonly_value


        venn3(subsets = (len(Lonly), len(Sonly), len(LS), len(Wonly), len(LW), len(SW), len(LSW)), set_labels = ('Large >3μm', 'Small 3-02μm', 'Whole water <0.22μm'), alpha = 0.5);

        plt.savefig("outputs/"+comm+"/D"+str(depths[d])+level+"_venn.png")
        plt.clf()
        plt.cla()
        plt.close()

    dfplot['Depth'] = dfplot['Depth'].astype(str)
    dfplot = dfplot.set_index('Depth')
    dfplot_normalized = dfplot/NewTotal *100

    return dfplot, dfplot_normalized, level


# In[ ]:


def calcperc_LSW(comm, separated, level):

    depths = [1, 5, 10, 30, 60]

    level = level

    dfplot = pd.DataFrame(columns=['Depth', 'NSF', 'LW', 'SW', 'LSW'])

    for d in range(len(depths)):
        sfd=separated[separated.depth==depths[d]]
        toptaxa = sfd[[level, 'feature_frequency', 'Taxon', 'size_code', 'weekn']].copy()

        toptaxa = toptaxa.drop_duplicates()
        df_agg = toptaxa.groupby(['size_code',level]).agg({'feature_frequency':sum})

        df_agg = df_agg.reset_index()
        resultpivot = df_agg.pivot_table(index=level, columns='size_code', values='feature_frequency')
        resultpivot = resultpivot.fillna(0)

        df = resultpivot.copy()

        Sonly = df[(df['L'] == 0) & (df['W'] == 0)]
        Wonly = df[(df['L'] == 0) & (df['S'] == 0)]
        Lonly = df[(df['S'] == 0) & (df['W'] == 0)]
        LW = df[(df['S'] == 0) & (df['W'] != 0) & (df['L'] != 0)]
        LS = df[(df['W'] == 0) & (df['S'] != 0) & (df['L'] != 0)]
        SW = df[(df['W'] != 0) & (df['S'] != 0) & (df['L'] == 0)]
        LSW = df[~(df == 0).any(axis=1)]

        total = df.to_numpy().sum()

        SFdf = Lonly, LS, Sonly
        SF = pd.concat(SFdf)
        SF_value = SF.to_numpy().sum()/total *100

        Sonly_value = Sonly.to_numpy().sum()/total *100


        Bothdf = LW, LSW, SW
        Both = pd.concat(Bothdf)
        Both_value = Both.to_numpy().sum()/total *100

        Wonly_value = Wonly.to_numpy().sum()/total *100

        Sonly_value = Sonly.to_numpy().sum()/total *100
        Lonly_value = Lonly.to_numpy().sum()/total *100
        LS_value = LS.to_numpy().sum()/total *100

        LW_value = LW.to_numpy().sum()/total *100
        SW_value = SW.to_numpy().sum()/total *100
        LSW_value = LSW.to_numpy().sum()/total *100

        NewTotal = Sonly_value + Lonly_value + LS_value + Wonly_value


        dfplot.loc[d,'Depth'] = depths[d]
        dfplot.loc[d,'NSF'] = Wonly_value
        dfplot.loc[d,'LW'] = LW_value
        dfplot.loc[d,'SW'] = SW_value
        dfplot.loc[d,'LSW'] = LSW_value


        venn3(subsets = (len(Lonly), len(Sonly), len(LS), len(Wonly), len(LW), len(SW), len(LSW)), set_labels = ('Large >3μm', 'Small 3-02μm', 'Whole water <0.22μm'), alpha = 0.5);

        plt.savefig("outputs/"+comm+"/D"+str(depths[d])+level+"_venn.png")
        plt.clf()
        plt.cla()
        plt.close()

    dfplot['Depth'] = dfplot['Depth'].astype(str)
    dfplot = dfplot.set_index('Depth')
    dfplot_normalized = dfplot/NewTotal *100

    return dfplot, dfplot_normalized, level


# In[ ]:


def calcperc_LS_W(comm, separated, level):

    depths = [1, 5, 10, 30, 60]

    level = level

    dfplot = pd.DataFrame(columns=['Depth', 'NSF', 'LW', 'SW'])

    for d in range(len(depths)):
        sfd=separated[separated.depth==depths[d]]
        toptaxa = sfd[[level, 'feature_frequency', 'Taxon', 'size_code', 'weekn']].copy()

        toptaxa = toptaxa.drop_duplicates()
        df_agg = toptaxa.groupby(['size_code',level]).agg({'feature_frequency':sum})

        df_agg = df_agg.reset_index()
        resultpivot = df_agg.pivot_table(index=level, columns='size_code', values='feature_frequency')
        resultpivot = resultpivot.fillna(0)

        df = resultpivot.copy()

        Sonly = df[(df['L'] == 0) & (df['W'] == 0)]
        Wonly = df[(df['L'] == 0) & (df['S'] == 0)]
        Lonly = df[(df['S'] == 0) & (df['W'] == 0)]
        LW = df[(df['S'] == 0) & (df['W'] != 0) & (df['L'] != 0)]
        LS = df[(df['W'] == 0) & (df['S'] != 0) & (df['L'] != 0)]
        SW = df[(df['W'] != 0) & (df['S'] != 0) & (df['L'] == 0)]
        LSW = df[~(df == 0).any(axis=1)]

        total = df.to_numpy().sum()

        SFdf = Lonly, LS, Sonly
        SF = pd.concat(SFdf)
        SF_value = SF.to_numpy().sum()/total *100

        Sonly_value = Sonly.to_numpy().sum()/total *100


        Bothdf = LW, LSW, SW
        Both = pd.concat(Bothdf)
        Both_value = Both.to_numpy().sum()/total *100

        Wonly_value = Wonly.to_numpy().sum()/total *100

        Sonly_value = Sonly.to_numpy().sum()/total *100
        Lonly_value = Lonly.to_numpy().sum()/total *100
        LS_value = LS.to_numpy().sum()/total *100

        LW_value = LW.to_numpy().sum()/total *100
        SW_value = SW.to_numpy().sum()/total *100
        LSW_value = LSW.to_numpy().sum()/total *100

        NewTotal = Sonly_value + Lonly_value + LS_value + Wonly_value


        dfplot.loc[d,'Depth'] = depths[d]
        dfplot.loc[d,'NSF'] = Wonly_value
        dfplot.loc[d,'LW'] = LW_value
        dfplot.loc[d,'SW'] = SW_value


        venn3(subsets = (len(Lonly), len(Sonly), len(LS), len(Wonly), len(LW), len(SW), len(LSW)), set_labels = ('Large >3μm', 'Small 3-02μm', 'Whole water <0.22μm'), alpha = 0.5);

        plt.savefig("outputs/"+comm+"/D"+str(depths[d])+level+"_venn.png")
        plt.clf()
        plt.cla()
        plt.close()

    dfplot['Depth'] = dfplot['Depth'].astype(str)
    dfplot = dfplot.set_index('Depth')
    dfplot_normalized = dfplot/NewTotal *100

    return dfplot, dfplot_normalized, level

def heatmap_top1(comm, sfd, level):
    if comm == 'chloroplast':
        level = 'PRSpecies'
    else:
        level = level

    toptaxa = sfd[['feature_id', 'feature_frequency', 'Taxon', 'size_code', 'depth', 'weekn', level]].copy()
    toptaxa.loc[toptaxa[level].isin(['Unassigned', 'uncultured+bacterium', 'g__uncultured']), level] = toptaxa['feature_id']
    toptaxa = toptaxa.drop_duplicates()

    df_agg = toptaxa.groupby(['size_code', level, 'depth', 'weekn']).agg({'feature_frequency': sum})
    topd = df_agg['feature_frequency'].groupby(['size_code', 'depth', 'weekn'], group_keys=False).nlargest(1).reset_index()

    topd[level] = topd[level].str.replace('g__', '', regex=False)

    genus_abundance = topd.groupby(level)['feature_frequency'].sum().sort_values(ascending=False)
    unique_genera = genus_abundance.index.tolist()
    unique_genera = sorted(unique_genera, reverse=True)

    max_colors = plt.get_cmap('tab20').N
    if len(unique_genera) > max_colors:
        top_genera = unique_genera[:max_colors - 1]
        topd[level] = topd[level].apply(lambda x: x if x in top_genera else 'Other')
        unique_genera = ['Other']+ sorted(top_genera, reverse=True)

    type_dic = {genus: i for i, genus in enumerate(unique_genera, start=1)}
    topd['comm_type'] = topd[level].map(type_dic)

    size_order = ['W', 'SL', 'S', 'L']
    depth_order = [f"{d}{s}" for d in ['1', '5', '10', '30', '60'] for s in size_order]
    topd["sc_weekn"] = topd["depth"].astype(str) + topd["size_code"]
    topd['sc_weekn'] = pd.Categorical(topd['sc_weekn'], categories=depth_order, ordered=True)

    topd = topd.sort_values(['sc_weekn'])
    glue = topd.pivot(index="sc_weekn", columns="weekn", values="comm_type")
    glue = glue[glue.columns].astype(float)

    cmap = plt.get_cmap('tab20', len(type_dic))

    sns.set_style('ticks')
    plt.figure(figsize=(5, 5))
    ax = sns.heatmap(glue, fmt='f', yticklabels=True, linewidths=.5, cmap=cmap)

    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([val + 0.5 for val in range(len(type_dic))])
    colorbar.set_ticklabels(unique_genera)

    ax.axhline(4, ls='--')
    ax.axhline(8, ls='--')
    ax.axhline(12, ls='--')
    ax.axhline(16, ls='--')

    ax.set_xticks(range(0, 16, 5))
    ax.set_ylabel("Depth and Size Code")
    ax.set_xlabel("Week Number")

    plt.savefig(f'outputs/{comm}/heatmap_top1_{level}.png', bbox_inches='tight', dpi=300)
    plt.show()




#subtitle = 'From Jan7 2022 to Apr27 2022'
def plot_stackedbar_p_SLNSF(comm, df, labels, colors, title, subtitle, level, xmax=110, xtick=10):

    fields = df.columns.tolist()

    #extract the df name as a string to label the output plot of this function
    tablelabel = f'{df=}'.split('=')[0]

    # figure and axis
    fig, ax = plt.subplots(1, figsize=(8, 5))
# plot bars
    left = len(df) * [0]
    for idx, name in enumerate(fields):
        plt.barh(df.index, df[name], left = left, color=colors[idx])
        left = left + df[name]
# title and subtitle
    plt.title(title, loc='left')
    #plt.text(0, ax.get_yticks()[-1] + 0.75, subtitle)
# legend
    plt.legend(labels, bbox_to_anchor=([1, 1, 0, 0]), ncol=1, frameon=False)
# remove spines
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
# format x ticks
    xticks = np.arange(0,xmax,xtick)
    xlabels = ['{}%'.format(i) for i in np.arange(0,xmax,xtick)]
    plt.xticks(xticks, xlabels)
# adjust limits and draw grid lines
    plt.ylim(-0.5, ax.get_yticks()[-1] + 0.5)
    ax.xaxis.grid(color='gray', linestyle='dashed')
    plt.gca().invert_yaxis()
    plt.ylabel("Depth (m)")
    plt.savefig('outputs/'+comm+'/'+level+'alldepths_stacked_perc_weighted'+ tablelabel +'.png', dpi=200, bbox_inches="tight")
    plt.show()

def uniqlist(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def grab_80perc(comm, df, perc, level, incl_mean=False):
    #count n ASVs that represent the top 80% of reads
    df.sort_values(by=['sampleid','ratio'], ascending=[True,False], inplace=True)
    df['cumul'] = df.groupby('sampleid')['ratio'].transform(pd.Series.cumsum)
    df_Grouped = df.groupby(['sampleid']).agg({'cumul': [('Count', 'count'),
                                                                   ('Topn80', lambda x: len(x[x<0.8]))]}).reset_index()

    #find the nearest larger neighbor to 80% (to include the asv that is within the 80%)
    df['difperc']=(df['cumul']-perc)
    df1=df[df['difperc'] > 0].groupby('sampleid').min().reset_index().reindex()

    #remove the ASVs not in the 80%
    df1 = df1[['sampleid', 'cumul']]
    df1.rename(columns={"cumul": "maxval"}, errors="raise", inplace=True)
    df = df.merge(df1,on='sampleid')
    df = df[df['cumul']<df['maxval']] #filter out the ASVs that are less than the 80% of reads

    #subselect only columns of interest
    df= df[['sampleid', level, 'feature_frequency', 'weekn', 'size_code', 'depth']]

    #get the list of top ASVs from the whole fraction
    dfW = df[df.size_code == 'W']
    df3 = dfW.groupby(['weekn', 'depth'])[level].apply(list).reset_index(name='tocompareagainst')
    df3['weekdpth'] = df3["weekn"].astype(str) + df3["depth"].astype(str)

    #get the list of top ASVs from all the size fractions
    dfSLSL = df[df.size_code != 'W']
    df4 = dfSLSL.groupby(['weekn', 'depth','size_code'])[level].apply(list).reset_index(name='tocomparewith')
    df4['tocomparewith'] = df4['tocomparewith'].apply(lambda x: list(set(x)))

    #merge with the list from Whole
    df5 = df4.merge(df3, on=['weekn', 'depth'], how='outer')

    #make a list of the top 80% reads ASVs ids
    df5['tocompareagainst'] = df5['tocompareagainst'].apply(lambda d: d if isinstance(d, list) else [])
    df5['tocompareagainst'] = df5['tocompareagainst'].apply(lambda x: list(set(x)))

    #make a list of the intersection, i.e. ASVs in top80% of reads shared between any size fraction and W
    df5['intersection'] = [list(set(a).intersection(set(b)))
                          for a, b in zip(df5.tocomparewith, df5.tocompareagainst)]

    #make a list of the differences
    df5['difference'] = df5['tocomparewith'].map(set) - df5['tocompareagainst'].map(set)

    #calculate the size of the intersection as a percetnage
    df5['lengthinter'] = df5['intersection'].str.len()
    df5['lengthtotal'] = df5['tocompareagainst'].str.len()
    df5['lengthperc'] = df5['lengthinter']/df5['lengthtotal']
    df5[["depth"]] = df5[["depth"]].apply(pd.to_numeric)
    df5.sort_values(['depth', 'size_code'], inplace=True)
    df5['weekdpth'] = df5["depth"].astype(str) + df5['size_code'].astype(str)

    #save the output as csv to inspect groups
    df5.to_csv('outputs/'+comm+'/intersection_with_W'+level+'_'+str(perc)+'.csv')

    ordered_yaxis = df5['weekdpth'].tolist()

    mylist = uniqlist(ordered_yaxis)


    glue = df5.pivot(index="weekdpth", columns="weekn", values="lengthperc")
    glue = glue.reindex(mylist)
    glue = glue[glue.columns].astype(float)
    if incl_mean == True:
        glue['mean'] = glue.mean(axis=1)
        glue.to_csv('outputs/'+comm+'/interesction_with_W_'+level+'_'+str(perc)+'.csv')

    sns.set_style('ticks')
    plt.figure(figsize=(5, 4.4))


    ax = sns.heatmap(glue, fmt='f', yticklabels=True, linewidths=.5, cmap=sns.color_palette("coolwarm", as_cmap=True))
    ax.axhline(3, ls='--')
    ax.axhline(6, ls='--')
    ax.axhline(9, ls='--')
    ax.axhline(12, ls='--')

    ax.set_xticks(range(0, 16, 5))

    plt.savefig('outputs/'+comm+'/interesction_with_W_'+level+'_'+str(perc)+'.png', bbox_inches='tight', dpi=300)

    plt.show()

def get_slopes(comm, df):
    forr2 = df[["nASVs", "weekn", 'sampleid', 'size_code', 'depth']].copy()
    forr2[["nASVs", "weekn"]] = forr2[["nASVs", "weekn"]].apply(pd.to_numeric)
    grouped_by_trims = forr2.groupby(['depth','size_code'])

    cols = forr2['size_code'].unique()

    linreg = (grouped_by_trims.apply(lambda x: pd.Series(stats.linregress(x['nASVs'], x['weekn'])))
            .rename(columns={
                0: 'slope',
                1: 'intercept',
                2: 'rvalue',
                3: 'pvalue',
                4: 'stderr'}))

    linreg['r2'] = linreg['rvalue']**2
    linreg = linreg.reset_index()


    tohm = linreg.pivot("depth", "size_code", "slope")
    tohm['mean'] = tohm.mean(axis=1)
    tohm['std'] = tohm.std(axis=1)
    tohm['CV'] = tohm['std']/tohm['mean']

    newcols=[]
    for col in cols:
        tohm['z'+col] = (tohm[col]-tohm['mean'])/tohm['std']
        newcols.append('z'+col)

    z_sc_df = tohm[newcols].copy()


    ax = sns.heatmap(tohm,cmap=("coolwarm"), annot=True, annot_kws={"fontsize":10}, cbar_kws={'label': 'Slope (ΔASVs/week'})
    plt.xlabel('Size fraction')
    plt.ylabel('Depth (m)')
    plt.savefig('outputs/'+comm+'/slopes.png', bbox_inches='tight', dpi=300)
    plt.clf()
    plt.cla()
    plt.close()

    ax = sns.heatmap(z_sc_df,cmap=("coolwarm"), annot=True, annot_kws={"fontsize":10}, cbar_kws={'label': 'Z scores'})
    plt.xlabel('Size fraction')
    plt.ylabel('Depth (m)')
    plt.savefig('outputs/'+comm+'/z_scores_slopes.png', bbox_inches='tight', dpi=300)
    plt.show()

    return tohm, z_sc_df



def run_RF(comm, depth, d_spc, newseparated):
    forrandomforest = d_spc[['feature_id', 'ratio', 'size_code','sampleid']]
    forrandomforest.drop_duplicates(inplace=True)
    RF = forrandomforest.pivot(index=['sampleid','size_code'],columns='feature_id', values='ratio').reset_index()
    RF.drop(columns='sampleid', inplace=True)
    RF = RF.fillna(0)

    RF['size_code'] = RF['size_code'].map({'S':1,'W':2, 'SL':3, 'L':4})

    # Split the data into features (X) and target (y)
    X = RF.drop('size_code', axis=1)
    y = RF['size_code']

    # Split the data into training and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

    rf = RandomForestClassifier()
    rf.fit(X_train, y_train)

    y_pred = rf.predict(X_test)

    accuracy = accuracy_score(y_test, y_pred)
    print("Accuracy:", accuracy)

    param_dist = {'n_estimators': randint(50,500),
                  'max_depth': randint(1,20)}

    # Create a random forest classifier
    rf = RandomForestClassifier()

    # Use random search to find the best hyperparameters
    rand_search = RandomizedSearchCV(rf, param_distributions = param_dist , n_iter=5, cv=5)

    # Fit the random search object to the data
    rand_search.fit(X_train, y_train)

    # Create a variable for the best model
    best_rf = rand_search.best_estimator_

    # Print the best hyperparameters
    print('Best hyperparameters:',  rand_search.best_params_)

    # Generate predictions with the best model
    y_pred = best_rf.predict(X_test)

    # Create the confusion matrix
    cm = confusion_matrix(y_test, y_pred)
    #ConfusionMatrixDisplay(confusion_matrix=cm).plot();

    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred, average='weighted')
    recall = recall_score(y_test, y_pred, average = 'weighted')

    print("Accuracy:", accuracy)
    print("Precision:", precision)
    print("Recall:", recall)

    # Create a series contain feature importances from the model and feature names from the training data
    feature_importances = pd.Series(best_rf.feature_importances_, index=X_train.columns).sort_values(ascending=False)

    feature_importances20 =feature_importances.head(20)
    #add taxonomic id to important features
    if comm == 'chloroplast':
        tax = newseparated[['feature_id','PRTaxon']].copy()
    else:
        tax = newseparated[['feature_id','Taxon', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']].copy()
    tax.drop_duplicates(inplace=True)
    new=feature_importances20.to_frame().reset_index()
    new = new.merge(tax, how='left', on='feature_id')
    new.to_csv('outputs/'+comm+'/top10predictors'+str(depth)+'.csv')
    return new, feature_importances

def notify():
    os.system('osascript -e \'display notification "Your code has finished running!" with title "Jupyter Notebook"\'')

#from chatgpt

def count_feature_id_presence_with_depth_and_W(directory_path, comm):
    """
    Counts feature ID presence, depth, and W sum for the given directory and community.

    Args:
    - directory_path: The base directory where the ANCOM files are located.
    - comm: The community type to locate the specific folder within ANCOM (e.g., chloroplast).

    Returns:
    - A list summarizing feature ID presence and W value information.
    """
    # Dictionary to track files, depths, and sum of W values for each feature_id
    feature_id_info = defaultdict(lambda: {"file_count": 0, "depths": set(), "files": set(), "W_sum": 0})

    # Regex pattern to match depth in file names and ensure "Trueonly" is present
    depth_pattern = re.compile(r"_D(1|5|10|30|60)_.*Trueonly")

    # Construct the path to the community folder
    search_path = os.path.join(directory_path, 'ANCOM', comm)
    #print(f"Searching in directory: {search_path}")

    # Walk through the community folder and all its subfolders to find all matching CSV files
    for root, dirs, files in os.walk(search_path):
        for file in files:
            if file.endswith(".csv"):
                #print(f"Considering file: {file}")  # Debugging: print all CSV files found
                if depth_pattern.search(file):
                    #print(f"File matches depth pattern: {file}")  # Debugging: print files that match depth pattern
                    file_path = os.path.join(root, file)

                    # Load CSV file
                    try:
                        df = pd.read_csv(file_path)
                    except Exception as e:
                        print(f"Error reading file {file_path}: {e}")
                        continue

                    # Check if 'feature_id' and 'W' columns exist
                    if 'feature_id' in df.columns and 'W' in df.columns:
                        #print(f"File {file_path} contains required columns")  # Debugging: Confirm file has required columns
                        # Get unique feature_ids in this file
                        for _, row in df.iterrows():
                            feature_id = row['feature_id']
                            W_value = row['W']

                            # Increment file count only once per file for each feature_id
                            if file_path not in feature_id_info[feature_id]["files"]:
                                feature_id_info[feature_id]["file_count"] += 1
                                feature_id_info[feature_id]["files"].add(file_path)

                            # Add depth to the set of depths for this feature_id
                            feature_id_info[feature_id]["depths"].add(f"Depth_{depth_pattern.search(file).group(1)}")

                            # Sum the W value
                            feature_id_info[feature_id]["W_sum"] += W_value
                    else:
                        print(f"Skipped {file_path}: 'feature_id' or 'W' column missing")

    # Convert to a list sorted by file count for easier readability
    feature_id_summary = sorted(
        [(feature_id, info["file_count"], sorted(info["depths"]), sorted(info["files"]), info["W_sum"])
         for feature_id, info in feature_id_info.items()],
        key=lambda x: x[1],
        reverse=True
    )

    #print(f"Collected feature_id_info: {dict(feature_id_info)}")  # Debugging: Print final collected data
    return feature_id_summary


def taxonomic_barplots(comm, table, depths, level, n=15, include_other=True):
    if comm == 'chloroplast':
        level = 'PRSpecies'

    # Filter data for the relevant depths
    sfd_all = table[table['depth'].isin(depths)]

    # Replace "Unassigned" or "uncultured+bacterium" at the given level with its 'feature_id'
    sfd_all[level] = sfd_all.apply(
        lambda row: row['feature_id'] if row[level] in ['Unassigned', 'uncultured+bacterium', 'c__uncultured', 'g__uncultured'] else row[level],
        axis=1
    )

    # Aggregate data globally to find the most abundant taxa across all depths
    global_abundance = sfd_all.groupby(level)['ratio'].sum().sort_values(ascending=False)

    # Select the top n most abundant taxa globally
    top_taxa = global_abundance.head(n - 1).index.tolist()

    if include_other:
        # Create a new column to group less abundant taxa into "Other"
        sfd_all['plot_taxa'] = sfd_all[level].apply(lambda x: x if x in top_taxa else 'Other')
        unique_taxa = ['Other'] + top_taxa
    else:
        sfd_all = sfd_all[sfd_all[level].isin(top_taxa)]
        sfd_all['plot_taxa']=sfd_all[level]
        unique_taxa = top_taxa

    # Format taxa labels for the legend (remove prefixes and "__")
    def format_taxa_label(taxon):
        return taxon.split("__")[-1]

    # Create a global discrete color palette for taxa
    n_taxa = len(unique_taxa)

    # Use a categorical colormap with up to 15 unique colors
    color_palette = plt.get_cmap('tab20').colors  # tab20 provides up to 20 colors
    color_dict = {taxon: color_palette[i % len(color_palette)] for i, taxon in enumerate(unique_taxa)}

    # Generate stacked bar plots for each depth
    for depth in depths:
        sfd = sfd_all[sfd_all['depth'] == depth]

        # Sort data by 'weekn' to ensure the correct order
        sfd = sfd.sort_values(by='weekn')

        size_codes = sfd['size_code'].unique()

        for idx, size_code in enumerate(size_codes):
            # Filter data for the current size_code
            top10d = sfd[sfd['size_code'] == size_code].copy()
            phyld = top10d.groupby(['weekn', 'date', 'plot_taxa'])['ratio'].sum().reset_index()

            all_weeks = pd.DataFrame({'weekn': range(1, 17)})
            phyld = phyld.merge(all_weeks, on='weekn', how='right')
            phyld.fillna({'ratio': 0}, inplace=True)
            phyld['date'].fillna('', inplace=True)  # or assign a placeholder value

            # Pivot the DataFrame to prepare for stacked bar plotting
            phyld_pivot = phyld.pivot(index=['weekn', 'date'], columns='plot_taxa', values='ratio').fillna(0).reset_index()

            # Determine x positions based on the 'weekn' order
            phyld_pivot = phyld_pivot.sort_values(by='weekn')  # Ensure data is ordered by weekn
            dates = phyld_pivot['date']  # Dates corresponding to the week numbers
            x = np.arange(len(dates))  # Positions for bars

            # Initialize stack heights
            bottom_stack = np.zeros(len(x))

            # Initialize the plot
            fig, ax1 = plt.subplots(figsize=(15, 7))

            # Iterate over each taxon and stack their relative abundances
            taxa_categories = phyld_pivot.columns[2:]  # Skip 'weekn' and 'date'
            for taxon in taxa_categories:
                # Aggregate taxon values by date
                heights = phyld_pivot[taxon].values

                # Plot the stacked bar for this taxon
                ax1.bar(
                    x,  # X positions
                    heights,  # Heights of the bars
                    bottom=bottom_stack,  # Start at the current stack height
                    alpha=0.8,
                    color=color_dict.get(taxon, 'gray')  # Use color from dictionary or default to gray
                )
                # Update the stack height
                bottom_stack += heights

            # Customize x-axis
            ax1.set_xticks(x[::2])
            ax1.set_xticklabels(dates[::2], fontsize=21) #rotation=45)
            ax1.set_ylabel("Relative Abundance", fontsize=24)
            ax1.set_xlabel("Date", fontsize=24)
            ax1.tick_params(axis='y', labelsize=21)

            # Do not add any legends to individual plots
            ax1.legend().remove()

            # Display the plot
            #plt.title(f"Stacked Bar Plot of Relative Abundance by Date for Size Code: {size_code}, Depth: {depth}")
            plt.tight_layout()
            plt.savefig(f'outputs/{comm}/stacked_bar_plot_size_code_{size_code}_depth_{depth}_{include_other}.png', dpi=300, bbox_inches='tight')
            plt.show()

    # Create a single legend figure for all depths
    fig_legend = plt.figure(figsize=(10, 1))
    legend_handles = [
        plt.Line2D([0], [0], color=color_dict[taxon], lw=10, label=format_taxa_label(taxon))
        for taxon in unique_taxa
    ]
    fig_legend.legend(
        handles=legend_handles,
        loc="center",
        ncol=1,
        title=level  # Use the level name as the legend title
    )
    plt.tight_layout()
    plt.savefig(f"outputs/{comm}/legend_for_all_depths.png", dpi=300, bbox_inches='tight')
    plt.show()


def plot_asv_depth_distribution(feature_id_summary):
    # Preparing data for plotting
    depths = ["Depth_1", "Depth_5", "Depth_10", "Depth_30", "Depth_60"]
    depth_counts = {depth: [] for depth in depths}
    feature_ids = []

    for entry in feature_id_summary:
        feature_id, _, entry_depths, _, _ = entry
        feature_ids.append(feature_id)

        # Assign depth presence (1 if present, 0 otherwise)
        for depth in depths:
            if depth in entry_depths:
                depth_counts[depth].append(1)
            else:
                depth_counts[depth].append(0)

    # Creating DataFrame for easy plotting
    df = pd.DataFrame(depth_counts, index=feature_ids)

    # Plotting a stacked bar chart
    ax = df.plot(kind="bar", stacked=True, figsize=(15, 8), cmap="viridis", edgecolor='none')

    # Adding labels and title
    plt.xlabel("Feature IDs (ASVs)")
    plt.ylabel("Depth Presence (1 = Present, 0 = Absent)")
    plt.title("Depth Distribution of ASVs")
    plt.xticks(ticks=range(len(feature_ids)), labels=feature_ids, rotation=90, fontsize=8)
    plt.legend(title="Depth", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()

# Function to filter top ASVs based on criteria
def filter_top_asvs(feature_id_summary, method="top_W_sum", n=50):
    if method == "top_W_sum":
        # Select top ASVs by highest W sum
        sorted_summary = sorted(feature_id_summary, key=lambda x: x[4], reverse=True)
        return sorted_summary[:n]
    elif method == "most_depths":
        # Select ASVs that appear in the most depths
        sorted_summary = sorted(feature_id_summary, key=lambda x: len(x[2]), reverse=True)
        return sorted_summary[:n]
    elif method == "random":
        # Randomly select ASVs
        return random.sample(feature_id_summary, min(n, len(feature_id_summary)))
    else:
        raise ValueError("Invalid method. Choose from 'top_W_sum', 'most_depths', or 'random'.")


#for plotting the metadata
def plot_nutrients(df, depth):
    if depth == 60:
        df = df.loc[df['depth'] == 60]
    else:
        df = df.loc[df['depth'] != 60]

    df = df[['weekn', 'depth', 'Phosphate', 'Silicate', 'Nitrate',
           'Ammonia', 'Chlorophyll A']]
    df.rename(columns={'Temperature_x': 'Temperature'},
              inplace= True)
    df = df.melt(id_vars=['weekn', 'depth'])

    nutrients = df[df['variable'] != 'Chlorophyll A']
    chlorophyll = df[df['variable'] == 'Chlorophyll A']

    maxvals = df.loc[df.groupby('variable')['value'].idxmax()]

    fig, ax1 = plt.subplots(figsize=(15, 7))
    sns.lineplot(
        data=nutrients,
        x="weekn", y="value",
        hue='variable', palette='viridis', ax=ax1
    )
    ax1.set_ylabel("Nutrients (Phosphate, Silicate, Nitrate, Ammonia)")
    ax1.axis(ymin=0,ymax=45)
    ax1.legend(title="Nutrients", bbox_to_anchor=(1.05, 1), loc='upper left')

    # Second y-axis for Chlorophyll A
    ax2 = ax1.twinx()
    sns.lineplot(
        data=chlorophyll,
        x="weekn", y="value",
        color='red', label='Chlorophyll A', ax=ax2
    )
    ax2.set_ylabel("Chlorophyll A")
    ax2.axis(ymin=0,ymax=14)
    ax2.legend(title="Chlorophyll A", bbox_to_anchor=(1.05, 0.8), loc='upper left')

    # Labels and save
    ax1.set_xlabel("Week Number")
    plt.title(f'Nutrients and Chlorophyll A in depths {depth}')
    plt.tight_layout()
    plt.savefig(f'outputs/{depth}_nutrients_chla.png')
    plt.show()

    return maxvals

# Updated heatmap function to show distinct W values for each depth considering only "Trueonly" files
def plot_asv_heatmap(comm, feature_id_summary, file_filter=None, directory_path=None, n=50):
    """
    Plots a heatmap of W values for ASVs across depths, considering only files that match the filter condition.

    Args:
    - comm: The community type, which is used to specify the directory for input files.
    - feature_id_summary: The feature summary generated by count_feature_id_presence_with_depth_and_W.
    - file_filter: Optional string to filter the filenames considered (e.g., "WSLSL" or "WSL").
                   If None, all files are considered.
    - directory_path: Optional base directory for locating the ANCOM files.

    Returns:
    - None (displays the heatmap)
    """
    # Filter top ASVs based on selected method
    top_asvs_summary = filter_top_asvs(feature_id_summary, method="top_W_sum", n=n)

    # Preparing data for heatmap
    depths = ["Depth_1", "Depth_5", "Depth_10", "Depth_30", "Depth_60"]
    taxonomic_assignments = {}
    w_values_dict = {}

    for entry in top_asvs_summary:
        feature_id, _, entry_depths, files, _ = entry
        taxonomic_label = feature_id  # Default to feature_id if no taxonomic assignment is found

        # Create dictionary to hold W values per depth for this ASV
        w_values_per_depth = {depth: 0 for depth in depths}

        for file_path in files:
            # Apply file filter condition if specified
            if file_filter and file_filter not in file_path:
                continue

            # Adjust the file path to include the community subdirectory
            if directory_path:
                full_file_path = os.path.join(directory_path, 'ANCOM', comm, os.path.basename(file_path))
            else:
                full_file_path = file_path

            # Extract depth from the file name
            try:
                depth_match = re.search(r"_D(1|5|10|30|60)_", full_file_path)
                if depth_match:
                    depth = f"Depth_{depth_match.group(1)}"
                else:
                    continue  # Skip files without depth information
            except Exception as e:
                print(f"Error extracting depth from file {full_file_path}: {e}")
                continue

            # Load the CSV to extract the correct W value and taxonomic assignment for that feature_id
            try:
                df = pd.read_csv(full_file_path)
                if 'feature_id' in df.columns and 'W' in df.columns:
                    # Sum the W values for the given feature_id in that file for the specific depth
                    w_value_sum = df[df['feature_id'] == feature_id]['W'].sum()

                    # Store W sum for the current depth
                    w_values_per_depth[depth] += w_value_sum

                    # Extract taxonomic assignment
                    if comm == 'chloroplast':
                        prtaxon = df[df['feature_id'] == feature_id]['PRTaxon'].values[0] if 'PRTaxon' in df.columns else None
                        if pd.notna(prtaxon):
                            levels = prtaxon.split('|')
                            last_identified = levels[-1]
                            if last_identified in ["Unassigned", "N/A", "", "unidentified"] and len(levels) > 1:
                                last_identified = levels[-2]
                            taxonomic_label = last_identified
                    else:
                        # Original logic for non-chloroplast data
                        genus = df[df['feature_id'] == feature_id]['Genus'].values[0] if 'Genus' in df.columns else None
                        species = df[df['feature_id'] == feature_id]['Species'].values[0] if 'Species' in df.columns else None
                        family = df[df['feature_id'] == feature_id]['Family'].values[0] if 'Family' in df.columns else None
                        taxon = df[df['feature_id'] == feature_id]['Taxon'].values[0] if 'Taxon' in df.columns else None

                        # Normalize values for comparison
                        def normalize_taxonomic_label(label):
                            return label if pd.notna(label) and label not in ["Unassigned", "N/A", "", "unidentified"] else None

                        genus = normalize_taxonomic_label(genus)
                        species = normalize_taxonomic_label(species)
                        family = normalize_taxonomic_label(family)

                        # Assign taxonomic label based on the given preference
                        if species and genus:
                            taxonomic_label = f"{genus} {species}"
                        elif genus:
                            taxonomic_label = genus
                        elif family:
                            genus_label = genus if genus else "unidentified"
                            taxonomic_label = f"{family} ({genus_label})"
                        elif taxon:
                            levels = taxon.split('; ')
                            last_identified = ""
                            next_unassigned = ""
                            for i, level in enumerate(levels):
                                if level.split("__")[1] not in ["Unassigned", "N/A", "", "unidentified"]:
                                    last_identified = level
                                else:
                                    next_unassigned = level
                                    break
                            if last_identified:
                                taxonomic_label = f"{last_identified}; {next_unassigned}" if next_unassigned else last_identified
            except Exception as e:
                print(f"Error reading file {full_file_path}: {e}")

        # Store the W values for each depth under the corresponding taxonomic label
        if taxonomic_label in w_values_dict:
            for depth in depths:
                w_values_dict[taxonomic_label][depth] += w_values_per_depth[depth]
        else:
            w_values_dict[taxonomic_label] = w_values_per_depth

    # Creating DataFrame for easy plotting
    heatmap_df = pd.DataFrame(w_values_dict).T  # Transpose to have taxonomic labels as rows
    heatmap_df.columns = depths

    # Plotting the heatmap
    plt.figure(figsize=(15, 10))
    sns.heatmap(heatmap_df, cmap="YlGnBu", annot=True, fmt=".1f", linewidths=0.5, cbar_kws={'label': 'Sum of W Values per Depth'})
    plt.xlabel("Depth")
    plt.ylabel("Taxonomic Assignments (Genus, Species, or Family)")
    plt.title("Heatmap of Top Taxonomic Sum of W Values Across Depths (Filtered by File Name)" if file_filter else "Heatmap of Top Taxonomic Sum of W Values Across Depths")
    plt.tight_layout()

    # Save the plot
    plt.savefig('outputs/'+comm+'/ancom_heatmap_plot.png', dpi=300)  # Save the plot to a file

    # Show the plot
    plt.show()

#this function was debugged by chatgpt:
    forMIMARKS= newseparated[["sampleid", "time_string_x", "size_code", "depth"]].drop_duplicates()
    size_fraction_mapping = {
        'L': '>3µm',
        'S': '0.2-3µm',
        'W': '>0.2µm',
        'P': 'Pooled',
        'SL': 'NA'
    }

    # Create the new column 'size_fraction' based on the mapping
def make_MIMARKS(newseparated):
    forMIMARKS= newseparated[["sampleid", "time_string_x", "size_code", "depth"]].drop_duplicates()
    size_fraction_mapping = {
        'L': '>3µm',
        'S': '0.2-3µm',
        'W': '>0.2µm',
        'P': 'Pooled',
        'SL': 'NA'
    }

    # Create the new column 'size_fraction' based on the mapping
    forMIMARKS = forMIMARKS[forMIMARKS['size_code'] != 'SL']
    forMIMARKS['size_fraction'] = forMIMARKS['size_code'].map(size_fraction_mapping)
    forMIMARKS.to_csv('forMIMARKS.csv', index=False)
    forMIMARKS = forMIMARKS[forMIMARKS['size_code'] != 'SL']
    forMIMARKS['size_fraction'] = forMIMARKS['size_code'].map(size_fraction_mapping)

    forMIMARKS.to_csv('forMIMARKS.csv', index=False)
