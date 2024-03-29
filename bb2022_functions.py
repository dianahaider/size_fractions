#!/usr/bin/env python
# coding: utf-8

# In[ ]:

#for importing, formatting and data manipulation
import pandas as pd
import numpy as np
import glob
import tempfile
from qiime2 import Artifact
import yaml
import json

#for plotting
import matplotlib, random
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from matplotlib.patches import Patch
import seaborn as sns
#sns.set(style="whitegrid")
import plotly.express as px
from IPython.display import display

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
import scipy.stats as stats
import statsmodels.api as sa
import statsmodels.formula.api as sfa
import scikit_posthocs as spPH

# Special thanks to Alex Manuele https://github.com/alexmanuele
def consolidate_tables(MG):
    if MG == '16S':
        comm = '02-PROKs'
    else :
        comm = '02-EUKs'

    table_list = glob.glob('{0}/table.qza'.format('/Users/Diana/Documents/escuela/phd/size_fractions/BB22_size-fraction-comparison-analysed/to_transfer/'+comm))
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

    #df.to_csv(comm+'/merged_all_tables.tsv', sep='\t', index=False)
    print("Successfully saved all tables.")
    return df, comm


# In[ ]:


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


# In[ ]:

def make_defract(all_md, separated):
    #make sure all size codes are indicated
    all_md["size_code"] = all_md["sampleid"].str.extract(r'[1-9][0-9]?[A-E]([L-S])')
    all_md["size_code"] = all_md["size_code"].fillna('W')

    #only keep values from weeks 1 to 16
    sep_SL = all_md[all_md.size_code != "W"]
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
    sep_SLRA = sepSLRA[separated.size_code != "W"]

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

    sep_S = separated[separated.size_code == "S"]
    sep_L = separated[separated.size_code == "L"]


    sep_WO.reset_index(inplace=True, drop=True)
    sep_SLRA.reset_index(inplace=True, drop=True)

    #newseparated = pd.concat([sep_SLRA.reset_index(drop=True), sep_WO.reset_index(drop=True)], axis=0).reset_index(drop=True)
    newseparated = pd.concat([sep_SLRA, sep_WO, sep_L, sep_S], ignore_index=True)

    newseparated['weekdepth'] = newseparated["weekn"].astype(str) + newseparated["depth"].astype(str)
    newseparated['avg'] = newseparated['nASVs'].groupby(newseparated['weekdepth']).transform('mean')
    newseparated['diff'] = newseparated['nASVs'] - newseparated['avg']

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


def pick_metadata(merged, depth='all', size_fraction='both', year='all', R='all', F='all', txsubset = 'all'):
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

    return separated


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

# In[ ]:


def taxbarplot(separated, level, depth, topn, palette_dict): #separated is the df, #level is a string of taxonomic level column name, depth is an integer
    sfd=separated[separated.depth==depth]
    toptaxa = sfd[['feature_id', 'feature_frequency', 'Taxon', 'size_code', 'depth','weekn', level]].copy()
    toptaxa = toptaxa.drop_duplicates()
    df_agg = toptaxa.groupby(['size_code',level, 'depth']).agg({'feature_frequency':sum})
    topd = df_agg['feature_frequency'].groupby('size_code', group_keys=False).nlargest(topn)
    topd = topd.to_frame()
    topd = topd.reset_index()


    df_agg = df_agg.reset_index()
    df_agg['set_name'] = df_agg['size_code']+df_agg['depth'].astype(str)

    cumulab = separated[['feature_frequency', 'depth', 'size_code', 'Genus']].copy()
    cumulab1 = cumulab.groupby(['Genus']).agg({'feature_frequency':sum})

    resultpivot = df_agg.pivot_table(index=level, columns='set_name', values='feature_frequency')
    resultpivot = resultpivot.fillna(0)
    resultpivot[resultpivot != 0] = 1
    tosave = pd.merge(resultpivot, cumulab1, left_index=True, right_index=True)
    tosave.to_csv(level+'_'+str(depth)+'16S_relab.csv')

    top10d_list = topd[level].unique()
    top10d = sfd.copy()
    top10d.loc[~top10d[level].isin(top10d_list), level] = 'Other' #isnot in top list
    phyld = top10d.groupby(['size_code','weekn', level])['ratio'].sum()
    phyld = phyld.reset_index()


    fig = px.bar(phyld, x="size_code", y="ratio", facet_col="weekn", color=level, labels={
                     "feature_frequency": "Relative abundance",
                     "size_code": "",
                     "weekn": "w"}, color_discrete_map=palette_dict)
    fig.update_xaxes(type='category', dtick=1)
    fig.update_layout(
        title="Relative abundance of top 10" + level + 'observed at Depth' + str(depth),
        yaxis_title="Relative abundance",
        xaxis_title="Size fraction",
        legend_title=level,
        font=dict(size=8)
    )

    fig.show()
    #fig.write_image("outputs/fig1.png")
    #fig.to_image(format="png")

    return phyld, top10d


# In[ ]:


def pcaplot(separated, depth, comm, columnperm, spc):

    if depth == 'all':
        df = separated.copy()
    else:
        df=separated[separated.depth==depth]


    if 'SL' in separated['size_code'].unique():
        #sizecode palette codes
        sizecodes = ['S', 'L', 'W', 'SL']
        palette_colors = sns.color_palette()
        palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}
        dicsc = pd.Series(df.size_code.values,index=df.sampleid).to_dict()
        color_rows_sc = {k: palette_dict[v] for k, v in dicsc.items()}
        seriescr = pd.Series(color_rows_sc)

    else:
        #sizecode palette codes
        sizecodes = ['S', 'L', 'W']
        palette_colors = sns.color_palette()
        palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}
        dicsc = pd.Series(df.size_code.values,index=df.sampleid).to_dict()
        color_rows_sc = {k: palette_dict[v] for k, v in dicsc.items()}
        seriescr = pd.Series(color_rows_sc)

    #month palette code
    df['Month'] = df['date'].str.split('-').str[1]
    months = ['Jan', 'Feb', 'Mar', 'May', 'Apr']
    palette_colors = sns.color_palette("flare")
    palette_dict_month = {monthname: color for monthname, color in zip(months, palette_colors)}
    dic = pd.Series(df.Month.values,index=df.sampleid).to_dict()
    color_rows_month = {k: palette_dict_month[v] for k, v in dic.items()}
    seriesmonthcr = pd.Series(color_rows_month)

    dfcolors = pd.DataFrame({'Month': seriesmonthcr,'Size code':seriescr})

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
        if weekNb < 10:
            return 'Pre-bloom'
        elif weekNb == 10 :
            return 'Bloom'
        elif weekNb > 10:
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
                       palette=palette_dict)#, size = 'Week_Group')#,palette=sns.color_palette("dark:salmon_r", as_cmap=True))
    plt.ylabel('PCo2 (' + str(pc2v) + '% variance explained)')
    plt.xlabel('PCo1 (' + str(pc1v) +'% variance explained)')
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
    for group, contributions in metadata_contributions.items():
        plt.barh(contributions, group) #range(1, len(contributions) + 1),

    plt.ylabel('Principal Component')
    plt.xlabel('Average Loading Contribution')
    sns.despine()
    plt.legend(frameon=False)
    plt.savefig('outputs/'+comm+'/D'+str(depth)+spc+'_PCAplot_brplot.png', dpi=200, bbox_inches="tight")
    plt.clf()
    plt.cla()
    plt.close()


    ##clustermap
    ax = sns.clustermap(distance_matrix, method="complete", cmap='RdBu', annot=True,
               yticklabels=True, row_colors = dfcolors,
               annot_kws={"size": 7}, figsize=(15,12));

    handles1 = [Patch(facecolor=palette_dict_month[key]) for key in palette_dict_month]
    plt.legend(handles1, palette_dict_month, title='Month',
               bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper left')

    plt.savefig('outputs/'+comm+'/D'+str(depth)+spc+'_clustermap.png', dpi=200, bbox_inches="tight")


    return pca, pca_features, sfdclr


# In[ ]:


def boxplot_depth(separated, comm, depth, ycolumn, yaxislabel='def'):

    if yaxislabel != 'def':
        ycol = ycolumn

    #sfd=separated[separated.depth==depth]
    sfd = separated.copy()

    #sfd_S = sfd[['size_code', 'nASVs', 'weekn']].copy()
    #sfd_S = sfd_S.drop_duplicates()
    #sdfpv = sfd_S.pivot(index='weekn', columns='size_code', values='nASVs')
    #fvalue, pvalue = stats.f_oneway(sdfpv['L'], sdfpv['S'], sdfpv['W'])

    sfd_LM = sfd[['size_code', 'nASVs']].copy()
    sfd_LM = sfd_LM.drop_duplicates()
    lm = sfa.ols('nASVs ~ C(size_code)', data=sfd_LM).fit()
    anova = sa.stats.anova_lm(lm)
    results = spPH.posthoc_ttest(sfd_LM, val_col='nASVs', group_col='size_code', p_adjust='holm')

    if 'SL' in separated['size_code'].unique():
        #sizecode palette codes
        sizecodes = ['S', 'L', 'W', 'SL']
        palette_colors = sns.color_palette()
        palette_dict = {sizecode: color for sizecode, color in zip(sizecodes, palette_colors)}

    else:
        #define color palettes
        sizecodes = ['S', 'L', 'W']
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


def upsetprep(comm, level, separated):

    depths = [1, 5, 10, 30, 60]

    cumulab = separated[['feature_frequency', 'depth', 'size_code', level]].copy()
    cumulab1 = cumulab.groupby([level]).agg({'feature_frequency':sum})

    for d in depths:
        #make csv
        sfd=separated[separated.depth==d]

        toptaxa = sfd[['feature_id', 'feature_frequency', 'Taxon', 'size_code', 'depth','weekn', level]].copy()
        toptaxa = toptaxa.drop_duplicates()
        df_agg = toptaxa.groupby(['size_code',level, 'depth']).agg({'feature_frequency':sum})
        topd = df_agg['feature_frequency'].groupby('size_code', group_keys=False).nlargest(10)
        topd = topd.to_frame()
        topd = topd.reset_index()

        df_agg = df_agg.reset_index()
        df_agg['set_name'] = df_agg['size_code']+df_agg['depth'].astype(str)

        resultpivot = df_agg.pivot_table(index=level, columns='set_name', values='feature_frequency')
        resultpivot = resultpivot.fillna(0)
        resultpivot[resultpivot != 0] = 1
        tosave = pd.merge(resultpivot, cumulab1, left_index=True, right_index=True)
        tosave.to_csv('csvs/'+comm+'/'+level+'_d'+str(d)+'_relab.csv')


        #make json
        data = {
            "file": "https://raw.githubusercontent.com/dianahaider/size_fractions/main/csvs/"+comm+'/'+level+'_d'+str(d)+'_relab.csv',
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

        with open('json/'+comm+'/'+level+'_d'+str(d)+'.json', 'w') as f:
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
    ax=sns.barplot(data=sfd_f, x="weekn", y="diff_p_id", hue="size_code", palette=palette_dict)#, hue="size_code")
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.title(ttl)
    plt.ylabel('Ratio difference')
    plt.xlabel('Week number')
    plt.savefig('outputs/'+comm+'/D'+str(depth)+fid+'.png', dpi=200, bbox_inches="tight")


# In[ ]:


def run_ancom(separated, sfdclr, depth, ancomcol):

    sfd=separated[separated.depth==depth]


    df_ancom = sfd[['sampleid', ancomcol]].copy()
    df_ancom = df_ancom.drop_duplicates()
    df_ancom = df_ancom.set_index('sampleid')

    results = ancom(table=sfdclr, grouping=df_ancom[ancomcol])

    DAresults = results[0].copy()
    DARejected_SC = DAresults.loc[DAresults['Reject null hypothesis'] == True]
    DARejected_SC.sort_values(by=['W'])

    taxonomy = sfd[['feature_id', 'Confidence', 'Taxon', 'Phylum', 'Class', 'Family', 'Genus', 'Species']].copy()
    taxonomy = taxonomy.drop_duplicates()
    DARejected_SC_taxonomy = pd.merge(DARejected_SC, taxonomy, on="feature_id", how="left")
    DARejected_SC_taxonomy.sort_values(by='W')

    prcentile = results[1].copy()

    return DARejected_SC_taxonomy, prcentile


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

    dfplot = pd.DataFrame(columns=['Depth', 'SF', 'NSF', 'Both'])

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

        dfplot.loc[d,'Depth'] = depths[d]
        dfplot.loc[d,'SF'] = SF_value
        dfplot.loc[d,'NSF'] = Wonly_value
        dfplot.loc[d,'Both'] = Both_value


        venn3(subsets = (len(Lonly), len(Sonly), len(LS), len(Wonly), len(LW), len(SW), len(LSW)), set_labels = ('Large >3μm', 'Small 3-02μm', 'Whole water <0.22μm'), alpha = 0.5);

        plt.savefig("outputs/"+comm+"/D"+str(depths[d])+level+"_venn.png")
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


# In[ ]:


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
