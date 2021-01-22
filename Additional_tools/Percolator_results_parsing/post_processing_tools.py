import os
import sys
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
import seaborn as sns
sns.set_style('whitegrid')
import numpy as np
import pandas as pd
import locale
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
locale.setlocale(locale.LC_ALL, '')

def parse_proteins_psm_results(input_file, output_file, output_columns=6):
    with open(input_file, 'r') as FR:
        with open(output_file, 'w') as FW:
            lines = FR.readlines()
            FW.write(lines[0])
            for line in lines[1:]:
                line = line.rstrip("\n")
                elements = line.split("\t")
                out = ("\t".join(elements[0:(output_columns-1)]))+"\t"+(",".join(elements[(output_columns-1):]))
                FW.write(out+"\n")
    
    return

def plot_psm_addition(can_psms,ext_psms,can_psms2,ext_psms2, label1='ms2ReScore search engine', label2='Unique MS2ReScore', label3='Prosit Andromeda only', label4='Prosit all features'):
        
    fig = plt.figure(figsize=[4,10])
    ax = fig.add_subplot(1,1,1)
    
    overlap = len(set(can_psms).intersection(set(ext_psms)))
    can_len = len(can_psms)-overlap
    ext_len = len(ext_psms)-overlap
    
    overlap2 = len(set(can_psms2).intersection(set(ext_psms2)))
    can_len2 = len(can_psms2)-overlap2
    ext_len2 = len(ext_psms2)-overlap2
    
    cans = [0,4]
    overlaps = [x+1 for x in cans]
    exts = [x+2 for x in cans]
    rects = ax.bar(cans, [can_len, can_len2])
    rects2 = ax.bar(overlaps, [overlap, overlap2])
    rects3 = ax.bar(exts, [ext_len, ext_len2])
    
    x=[0,1,2,4,5,6]
    y=[can_len,overlap,ext_len,can_len2,overlap2,ext_len2]
    for i in range(6):
        ax.text(x[i]-0.1, y[i]+12000, "{:,d}".format(y[i]).replace(","," "),rotation=90)
    
    plt.xticks([0,1,2,4,5,6], [label1, 'Overlap', label2,
                              label3, 'Overlap', label4],
               rotation=90)
    plt.grid(axis='x')
    plt.tight_layout()
    sns.despine(right=True, left=True, bottom=False)
    
    plt.show()
    
    return

def plot_overlap_ms2rescore_prosit(ms2rescore_psms, prosit_psms, label1='MS2ReScore', label2='Prosit'):
    
    fig = plt.figure(figsize=[2,10])
    ax = fig.add_subplot(1,1,1)
    
    overlap = len(set(ms2rescore_psms).intersection(set(prosit_psms)))
    ms2rescore_len = len(ms2rescore_psms)-overlap
    prosit_len = len(prosit_psms)-overlap
    
    rects = ax.bar(0, ms2rescore_len)
    rects2 = ax.bar(1, overlap)
    rects3 = ax.bar(2, prosit_len)
    
    plt.xticks([0,1,2], [label1, 'Overlap',label2],
               rotation=90)
    
    plt.grid(axis='x')
    plt.tight_layout()
    sns.despine(right=True, left=True, bottom=False)
    
    plt.show()
    
    return

def split_ms2rescore_psmid(df, id_column="PSMId"):

    df['Raw file'], df['Scan number'] = zip(*df.apply(ms2rescore_id_splitter, id_column=id_column, axis=1))
    
    return df

def ms2rescore_id_splitter(row, id_column="PSMId"):

    els = row[id_column].split(".")
    
    return(els[0], els[1])

def split_prosit_psmid(df, id_column="PSMId"):
    
    df['Raw file'], df['Scan number'], df['Charge'], df['Scan event'] = zip(*df.apply(prosit_id_splitter, id_column=id_column, axis=1))
    
    return df

def prosit_id_splitter(row, id_column="PSMId"):
    
    els = row[id_column].split("-")
    
    return(els[0]+"-"+els[1], els[2], els[4], els[5])

def ms2rescore_make_commonid(row, pepcolumn='ModPeptide'):
    
    common_id=row['Raw file']+"-"+str(row["Scan number"])+"-"+row[pepcolumn].replace("(ox)","").replace("(ac)","").replace("_","")+"-"+str(row["ChargeN"])
    
    return common_id

def prosit_make_commonid_df(df):
    
    df['CommonId'] = df.apply(prosit_make_commonid, axis=1)
    
    return df

def prosit_make_commonid(row, pepcolumn="Peptide"):
    
    common_id = row['Raw file']+"-"+str(row['Scan number'])+"-"+row[pepcolumn].replace("_","").replace(".","").replace("(ac)","").replace("(ox)","").replace("m","M")+"-"+str(row["Charge"])
    
    return common_id

def msmstxt_make_commonid(df):
    
    df['CommonId'] = df.apply(msmstxt_make_commonid_row, axis=1)
    
    return df

def msmstxt_make_commonid_row(row):
    
    common_id = row['Raw file']+"-"+str(row['Scan number'])+"-"+row['Sequence'].replace("(ac)","").replace("(ox)","")+"-"+str(row['Charge'])
    
    return common_id

def prosit_pepseq_converter_row(row):
    
    pepseq = row['Peptide'].replace("_","").replace(".","").replace("m","M")
    modpepseq = row['Peptide'].replace("_","").replace(".","").replace("m","M(ox)")
    lenseq = len(pepseq)
    
    return(pepseq, lenseq, modpepseq)

def prosit_pepseq_converter(df):
    
    df['PepSeq'], df['SeqLen'], df['ModPepSeq'] = zip(*df.apply(prosit_pepseq_converter_row, axis=1))
    
    return df

def ms2rescore_pepseq_converter_row(row):
    
    pepseq = row['ModPeptide'].replace("_","").replace(".","").replace("(ox)","").replace("(ac)","")
    lenseq = len(pepseq)
    
    return(pepseq, lenseq)

def ms2rescore_pepseq_converter(df):
    
    df['PepSeq'], df['SeqLen'] = zip(*df.apply(ms2rescore_pepseq_converter_row, axis=1))
    
    return df

def plot_weights(weights_file, title="weights", figure_width=16):
    
    with plt.style.context('seaborn-whitegrid'):
        weights = pd.read_csv(weights_file, sep="\t").iloc[::3]
        weights.iloc[:,:] = weights.iloc[:,:].astype('float64')
        weights = weights.apply(np.mean, axis=0).sort_values()
        abs_weights = abs(weights)
        ext_weights = max(abs_weights)
        colorbar_weights = [abs_weights, -abs_weights]
        colorbar_weights = pd.concat(colorbar_weights)

        fig = plt.figure(figsize=[figure_width,10])
        ax = fig.add_subplot(1,1,1)

        colors = cm.RdBu((weights+max(abs_weights)) / float(ext_weights*2))
        plot = plt.scatter(colorbar_weights, colorbar_weights, c = colorbar_weights, cmap = 'RdBu')
        plt.tick_params(labelsize=11)
        plt.clf()
        plt.colorbar(plot)
        plt.ylabel('Weight')

        #plt.title(title)

        weights.plot.bar(color=colors)
        
        plt.savefig("weights.svg", format="svg")
    
    return

def write_purified_baseline_pin(input_pin, outputfile='output.pin',features=['RawScore','RawDeltaScore']):
    
    #Define output columns
    output_pre_cols=['SpecId','Label','ScanNr']
    output_post_cols=['ChargeN','enzInt','ModPepti, size=18de','Proteins']
    output_cols = output_pre_cols+features+output_post_cols
    
    #Restructure
    output_df = input_pin.loc[:,output_cols]
    output_df.loc[:,'Proteins'] = output_df.loc[:,'Proteins'].apply(lambda x: x.replace(","," "))
    
    #Write output
    with open(outputfile, 'w') as FW:
        output_df.to_csv(FW, sep="\t", header=True, index=False)
    
    return output_df

def fdr_plot(dataframes, q_value_column='q-value', title='FDR plot', underlimit_x=0.00001):
    
    with plt.style.context('seaborn-darkgrid'):
    
        #Init figure and axes
        fig = plt.figure(figsize=[10,6])
        ax = fig.add_subplot(1,1,1)

        for label, df in dataframes.items():
            #Sort based on q values
            df_sorted_qval = df.sort_values(q_value_column, axis=0).copy().reset_index(drop=True)
            #Drop decoys
            df_sorted_qval['cum_count'] = (df_sorted_qval['is_decoy']==False).cumsum()
            #Make plot
            ax.plot(df_sorted_qval[q_value_column], df_sorted_qval['cum_count'], label=label)
            if label=='Andromeda + Percolator':
                plt.axhline(y=min(df_sorted_qval[df_sorted_qval[q_value_column]>=0.01].loc[:,'cum_count']), color='grey', linestyle='--', linewidth=0.75)
                plt.axvline(x=0.01,color='grey', linestyle='--', linewidth=0.75)
                ax.text(0.0000011, 170000, 'Andromeda 1% FDR threshold', fontsize=16)
        
        #Figure details
        plt.xlabel('FDR (log)', size=18)
        plt.ylabel('Number of identified spectra', size=18)
        plt.xscale("log")
        plt.xlim(underlimit_x,1)
        plt.tick_params(labelsize=14)
        plt.legend(fontsize=16, loc='lower right')
        #plt.title(title, size=20)
        plt.savefig("fdr_plot_prosit_ms2rescore.svg", format="svg")
        plt.show()
    
    return

def fdr_plot_expol(dataframes, q_value_column='q-value-expol', title='FDR plot', underlimit_x=0.00001):
    
    #Init figure and axes
    fig = plt.figure(figsize=[10,6])
    ax = fig.add_subplot(1,1,1)
    
    for label, df in dataframes.items():
        #Sort based on q values
        df_sorted_qval = df.sort_values(q_value_column, axis=0).copy().reset_index(drop=True)
        #Drop decoys
        df_sorted_qval['cum_count'] = (df_sorted_qval['is_decoy']==False).cumsum()
        #Make plot
        ax.plot(df_sorted_qval[q_value_column], df_sorted_qval['cum_count'], label=label)
    
    #Figure details
    plt.xlabel('FDR (log scale)')
    plt.ylabel('Number of identified spectra')
    plt.xscale("log")
    plt.xlim(underlimit_x,1)
    plt.tick_params(labelsize=14)
    plt.legend(fontsize=16)
    plt.title(title, size=20)
    plt.show()
    
    return

def tp_plot(dataframes, q_value_column='q-value', title='True positives plot', underlimit_x=0.00001):
    
    with plt.style.context('seaborn-darkgrid'):
        #Init figure and axes
        fig = plt.figure(figsize=[10,6])
        ax = fig.add_subplot(1,1,1)

        for label, df in dataframes.items():
            #Sort based on q values
            df_sorted_qval = df.sort_values(q_value_column, axis=0).copy().reset_index(drop=True)
            #Drop decoys
            df_sorted_qval['cum_count'] = (df_sorted_qval['is_decoy']==False).cumsum()
            df_sorted_qval['tp_count'] = df_sorted_qval['cum_count'] - df_sorted_qval[q_value_column]*df_sorted_qval['cum_count']
            #Make plot
            ax.plot(df_sorted_qval[q_value_column], df_sorted_qval['tp_count'], label=label)
            if label=='Andromeda + Percolator':
                plt.axhline(y=min(df_sorted_qval[df_sorted_qval[q_value_column]>=0.01].loc[:,'tp_count']), color='grey', linestyle='--', linewidth=0.75)
                plt.axvline(x=0.01,color='grey', linestyle='--', linewidth=0.75)
                ax.text(0.0000011, 170000, 'Andromeda 1% FDR threshold', fontsize=16)
               

        #Figure details
        plt.xlabel('FDR (log)', size=18)
        plt.ylabel('Number of true positives', size=18)
        plt.xscale("log")
        plt.xlim(underlimit_x,1)
        plt.tick_params(labelsize=14)
        plt.legend(fontsize=16)
        #plt.title(title, size=20)
        plt.savefig("tp_plot_prosit_ms2rescore.svg", format="svg")
        plt.show()
    
    return

def extrapolate_qvals(df, start_train_x=10**(-2)):
    
    #Get training data
    train_data = df[(df['q-value']< start_train_x) & (df['q-value']>(10**(-4)))]
    train_x = np.array(train_data.score).reshape((-1, 1))
    train_y = np.log10(np.array(train_data.loc[:,'q-value']))
    
    #Model
    model = LinearRegression()
    model.fit(train_x, train_y)
    r_sq = model.score(train_x, train_y)
    print('coefficient of determination:', r_sq)
    print('intercept:', model.intercept_)
    print('slope:', model.coef_)
    
    #Resubstitute in df
    df['q-value-expol'] = df['q-value'].copy()
    
    def change_qval(row, model=model):
        if(row['q-value']<(10**(-4))):
            q_val_expol = 10**float(model.predict(np.array(row['score']).reshape((-1,1))))
        else:
            q_val_expol = row['q-value']
        return q_val_expol
    
    df['q-value-expol'] =  df.apply(change_qval, model=model, axis=1)
    
    return df

def extrapolate_qvals_poly(df, start_train_x=10**(-2)):
    
    #Get training data
    train_data = df[(df['q-value']< start_train_x) & (df['q-value']>(10**(-4)))]
    train_x = np.array(train_data.score).reshape((-1, 1))
    train_y = np.log10(np.array(train_data.loc[:,'q-value']))
    
    #Transformer
    transformer = PolynomialFeatures(degree=2, include_bias=False)
    transformer.fit(train_x)
    train_x_ = transformer.transform(train_x)
    
    #Model
    model = LinearRegression()
    model.fit(train_x_, train_y)
    r_sq = model.score(train_x_, train_y)
    print('coefficient of determination:', r_sq)
    print('intercept:', model.intercept_)
    print('slope:', model.coef_)
    
    #Resubstitute in df
    df['q-value-expol'] = df['q-value'].copy()
    
    def change_qval(row, model=model, transformer=transformer):
        if(row['q-value']<(10**(-4))):
            pred_x = np.array(row['score']).reshape((-1,1))
            pred_x_ = transformer.transform(pred_x)
            q_val_expol = 10**float(model.predict(pred_x_))
        else:
            q_val_expol = row['q-value']
        return q_val_expol
    
    df['q-value-expol'] =  df.apply(change_qval, model=model, transformer=transformer, axis=1)
    
    return df


def plot_percolatorscore_qval(dataframes, q_value_column='q-value', percolator_score_column='score', xlim="None", log_x=False, ylim="None"):
    
    fig = plt.figure(figsize=[10,6])
    ax = fig.add_subplot(1,1,1)
    
    for label, df in dataframes.items():
        df_sorted_score = df.sort_values(percolator_score_column, axis=0).copy().reset_index(drop=True)
        
        if ylim!="None":
            df_sorted_score = df_sorted_score[(df_sorted_score[q_value_column]>ylim[0]) & (df_sorted_score[q_value_column]<ylim[1])]
        ax.plot(df_sorted_score[percolator_score_column],df_sorted_score[q_value_column], label=label)
    
    plt.ylabel('q-value (log scale)')
    if log_x:
        plt.xlabel('Percolator score (log)')
    else:
        plt.xlabel('Percolator score')
    if log_x:
        plt.xscale("log")
    plt.yscale("log")
    if xlim!="None":
        plt.xlim(xlim)
    plt.legend()
    plt.show()
    
    return

def feature_target_decoy_hist(pin, feature_name='spectral_angle', title='Density histogram', xlabel='', binsize=10, zoom=False, zoomval=2000):
    
    with plt.style.context('seaborn-darkgrid'):
        #Init figure
        fig = plt.figure(figsize=[10,6])
        ax = fig.add_subplot(1,1,1)

        targets = pin[pin['Label']==1].loc[:,feature_name]
        decoys = pin[pin['Label']==-1].loc[:,feature_name]

        plt.hist(targets, label='Targets', color='lightblue', bins=np.arange(0,int(math.ceil(max(targets))),binsize))
        plt.hist(decoys, label='Decoys', color='orange', bins=np.arange(0,int(math.ceil(max(decoys))),binsize))
        #sns.distplot(targets, label='Targets', hist=True, kde=False, bins=50)
        #sns.distplot(decoys, label='Decoys', hist=True, kde=False, bins=50)

        if zoom:
            bottom, top = plt.ylim()
            plt.ylim(bottom, zoomval)

        plt.ylabel('Number of PSMs', size=18)
        plt.xlabel(xlabel, size=18)
        plt.tick_params(labelsize=14)
        plt.legend()
        plt.savefig("feature_target_decoy_hist.svg", format="svg")
        plt.legend(fontsize=16)
        plt.show()
        
    return

def post_score_target_decoy_hist(perc_results, title='Density histogram', binsize=1, score_column='score', y2='auto'):
    
    with plt.style.context('seaborn-darkgrid'):
        #Init figure
        fig = plt.figure(figsize=[10,6])
        ax = fig.add_subplot(1,1,1)

        targets = perc_results[perc_results['is_decoy']==False].loc[:,score_column]
        decoys = perc_results[perc_results['is_decoy']==True].loc[:,score_column]
        xlabel=''

        if score_column=='score':
            plt.hist(targets, label='Targets', color='lightblue', bins=np.arange(-5, int(math.ceil(max(targets))), binsize))
            plt.hist(decoys, label='Decoys', color='orange', bins=np.arange(-5, int(math.ceil(max(decoys))), binsize))
            xlabel='Percolator score'
        elif score_column=='posterior_error_prob':
            plt.hist(targets, label='Targets', color='lightblue', bins=10**np.arange(-15,0,binsize))
            plt.hist(decoys, label='Decoys', color='orange', bins=10**np.arange(-15,0,binsize))
            #plt.yscale("log")
            plt.xscale("log")
            xlabel='Posterior error probability'
        if score_column=='score':
            plt.xlim(-5,10)
        plt.ylabel('Number of PSMs', size=18)
        plt.xlabel(xlabel, size=18)
        if y2!='auto':
            plt.ylim((0, float(y2)))
        
        #plt.title(title)
        plt.legend()
        plt.tick_params(labelsize=14)
        plt.legend(fontsize=16)
        if score_column=='score':
            plt.savefig("percolator_hist.svg", format="svg")
        elif score_column=='posterior_error_prob':
            plt.savefig("pep_hist.svg", format="svg")

        plt.show()
    
    return

def feature_target_decoy_dist(pin, feature_name='spectral_angle', title='Density distribution'):
    
    #Init figure
    fig = plt.figure(figsize=[10,6])
    ax = fig.add_subplot(1,1,1)
    
    targets = pin[pin['Label']==1].loc[:,feature_name]
    decoys = pin[pin['Label']==-1].loc[:,feature_name]
    
    sns.distplot(targets, label='Targets', hist=False, kde=True)
    sns.distplot(decoys, label='Decoys', hist=False, kde=True)
    
    plt.ylabel('Density')
    plt.title(title)
    plt.legend()
    
    plt.show()
    
    return
    
def features_joint_plot(pin, feature1='andromeda', feature2='spectral_angle'):
    
    targets = pin[pin['Label']==1]
    decoys = pin[pin['Label']==-1]
    
    make_joint_plot(targets, decoys, feature1=feature1, feature2 = feature2)
    
    return

def mixed_features_joint_plot(pin1, pin2, feature1='spec_pearson_norm', feature2='spectral_angle'):
    
    merged_pin = pd.merge(pin1, pin2, how='inner', on=['CommonId', 'Raw file', 'Scan number','Label'])
    targets = merged_pin[merged_pin['Label']==1]
    decoys = merged_pin[merged_pin['Label']==-1]
    
    make_joint_plot(targets, decoys, feature1=feature1, feature2=feature2)
    
    return

def make_joint_plot(targets, decoys, feature1='andromeda', feature2='spectral_angle'):
    
    sns.set_style('darkgrid')
    sns.set(font_scale=1.4)
    graph = sns.JointGrid(x=feature1, y=feature2, data=targets)
    
    graph.x = targets.loc[:,feature1]
    graph.y = targets.loc[:,feature2]
    graph.plot_joint(plt.scatter, marker='.', c='orange', alpha=0.2, s=5)
    #Legend for scatter
    orange_patch = mpatches.Patch(color='orange', label='Targets')
    red_patch = mpatches.Patch(color='red', label='Decoys')
    plt.legend(handles=[orange_patch, red_patch], fontsize=14)
    graph.plot_marginals(sns.kdeplot, color="orange", shade=True, shade_lowest=True, legend=False)

    graph.x = decoys.loc[:,feature1]
    graph.y = decoys.loc[:,feature2]
    handle_decoy = graph.plot_joint(plt.scatter, marker='x', c='red',alpha=0.2, s=5)
    graph.plot_marginals(sns.kdeplot, color="red", shade=True, shade_lowest=True, legend=False)
    
    
    plt.gcf().set_size_inches(10, 10)
    for ax in plt.gcf().axes:
        l = ax.get_xlabel()
        ax.set_xlabel(l, fontsize=18)
        l = ax.get_ylabel()
        ax.set_ylabel(l, fontsize=18)
    if feature1=='spec_pearson_norm':
        graph.ax_joint.set_xlabel('Pearson correlation (MS2ReScore)')
    elif feature1=='cos':
        graph.ax_joint.set_xlabel('Cosine score (MS2ReScore)')
    if feature2=='spectral_angle':
        graph.ax_joint.set_ylabel('Spectral angle (Prosit)')
    if feature2=='spec_pearson_norm':
        graph.ax_joint.set_ylabel('Pearson correlation (MS2ReScore)')
    if feature2=='cos':
        graph.ax_joint.set_ylabel('Cosine score (MS2ReScore)')
    if feature1=='andromeda' or feature1=='RawScore':
        graph.ax_joint.set_xlabel('Andromeda score')
    sns.despine()
    plt.tight_layout()
    plt.savefig('joint_plot.png', format="png", dpi=300)
    plt.savefig('joint_plot.tiff', format="tiff", dpi=300)
    plt.show()
    
    return

def spectrum_mirror_plot(ms2rescore_spectrum_info):
    
    with plt.style.context('seaborn-darkgrid'):
        fig, axes = plt.subplots(figsize=(15, 5))

        axes.vlines(x=ms2rescore_spectrum_info['mz'], ymin=-0.0001, ymax=-ms2rescore_spectrum_info['target'], color='orange', label='Empirical spectrum')
        axes.vlines(x=ms2rescore_spectrum_info['mz'], ymin=0.0001, ymax=ms2rescore_spectrum_info['prediction'], color='purple', label='MS2PIP prediction')

        #Set 
        axes.set_xlabel('m/z', size=18)
        axes.set_ylabel('Intensity', size=18)

        #Make absolute y ticks
        ticks =  axes.get_yticks()
        ticks = axes.set_yticklabels(['{:.0f}'.format(abs(tick)) for tick in ticks])
        plt.legend(loc='center left', bbox_to_anchor=(0.8, 1.1))
        plt.tight_layout()
        plt.savefig("mirror_plot.svg", format="svg")
        plt.show()
    
    return


def spectrum_difference_plot(ms2rescore_spectrum_info):
    
    with plt.style.context('seaborn-darkgrid'):
        ms2rescore_spectrum_info.target = ms2rescore_spectrum_info.target.abs()
        ms2rescore_spectrum_info.prediction = ms2rescore_spectrum_info.prediction.abs()

        yions_common_parts = pd.DataFrame(index=ms2rescore_spectrum_info.index, columns=['m/z', 'intensity'])
        bions_common_parts = pd.DataFrame(index=ms2rescore_spectrum_info.index, columns=['m/z', 'intensity'])
        experimental_is_bigger = pd.DataFrame(index=ms2rescore_spectrum_info.index, columns=['m/z', 'intensity'])
        prediction_is_bigger = pd.DataFrame(index=ms2rescore_spectrum_info.index, columns=['m/z', 'intensity'])
        for i in ms2rescore_spectrum_info.index:
            if ms2rescore_spectrum_info.loc[i,'target']<ms2rescore_spectrum_info.loc[i,'prediction']:
                if ms2rescore_spectrum_info.loc[i,'ion']=='B':
                    bions_common_parts.loc[i,'m/z'] = ms2rescore_spectrum_info.loc[i,'mz']
                    bions_common_parts.loc[i,'intensity'] = ms2rescore_spectrum_info.loc[i,'target']
                else:
                    yions_common_parts.loc[i,'m/z'] = ms2rescore_spectrum_info.loc[i,'mz']
                    yions_common_parts.loc[i,'intensity'] = ms2rescore_spectrum_info.loc[i,'target']
                prediction_is_bigger.loc[i,'m/z'] = ms2rescore_spectrum_info.loc[i,'mz']
                prediction_is_bigger.loc[i,'intensity'] = ms2rescore_spectrum_info.loc[i,'prediction']
            else:
                if ms2rescore_spectrum_info.loc[i,'ion']=='B':
                    bions_common_parts.loc[i,'m/z'] = ms2rescore_spectrum_info.loc[i,'mz']
                    bions_common_parts.loc[i,'intensity'] = ms2rescore_spectrum_info.loc[i,'prediction']
                else:
                    yions_common_parts.loc[i,'m/z'] = ms2rescore_spectrum_info.loc[i,'mz']
                    yions_common_parts.loc[i,'intensity'] = ms2rescore_spectrum_info.loc[i,'prediction']
                experimental_is_bigger.loc[i,'m/z'] = ms2rescore_spectrum_info.loc[i,'mz']
                experimental_is_bigger.loc[i,'intensity'] = ms2rescore_spectrum_info.loc[i,'target']

        fig, axes = plt.subplots(figsize=(15, 8))
        axes.vlines(x=experimental_is_bigger['m/z'], ymin=-0.0001, ymax=experimental_is_bigger['intensity'], color='orange', label='Empirical spectrum excess', linewidth=2)
        axes.vlines(x=prediction_is_bigger['m/z'], ymin=-0.0001, ymax=prediction_is_bigger['intensity'], color='purple', label='Prediction excess', linewidth=2)
        axes.vlines(x=bions_common_parts['m/z'], ymin=0.0001, ymax=bions_common_parts['intensity'], color='lightblue', label='b-ions', linewidth=2)
        axes.vlines(x=yions_common_parts['m/z'], ymin=0.0001, ymax=yions_common_parts['intensity'], color='red', label='y-ions', linewidth=2)

        #Set 
        axes.set_xlabel('m/z', size=18)
        axes.set_ylabel('Intensity', size=18)
        #Make absolute y ticks
        ticks =  axes.get_yticks()
        ticks = axes.set_yticklabels(['{:.0f}'.format(abs(tick)) for tick in ticks])
        plt.legend(loc='center left', bbox_to_anchor=(0.75, 1.1))
        plt.tight_layout()
        plt.savefig("spectrum_difference_plot.svg", format="svg")
        plt.show()

    return

def gain_loss_plot(prosit_and, prosit_combi, tech='MS$^{2}$ReScore'):
    
    #Init
    fdr_ths = np.logspace(-1, -5, num=17)
    and_fdr_th = 0.01
    overlaps=[]
    losts=[]
    gains=[]
    overlap_at_001=''
    
    #Andromeda reference situation
    and_common_ids = prosit_and[(prosit_and['is_decoy']==False) & (prosit_and['q-value']<=and_fdr_th)].CommonId
    and_len_ids=len(and_common_ids)
    
    #Calculate
    for fdr_th in fdr_ths:
        prosit_common_ids = prosit_combi[(prosit_combi['is_decoy']==False) & (prosit_combi['q-value']<=fdr_th)].CommonId
        overlap = len(set(and_common_ids).intersection(set(prosit_common_ids)))/and_len_ids
        overlaps.append(overlap)
        losts.append(1-overlap)
        gain = len(set(prosit_common_ids) - set(and_common_ids))/and_len_ids
        gains.append(gain)
        if fdr_th==0.01:
            overlap_at_001=overlap
    
    #Plot
    with plt.style.context('seaborn-darkgrid'): 
        #Init fig
        fig = plt.figure(figsize=[10,6])
        ax = fig.add_subplot(1, 1, 1)
        
        #Plot
        ax.plot(fdr_ths, overlaps, '-o', color='orange', label='In common with Andromeda +\n Percolator at $10^{-2}$')
        ax.plot(fdr_ths, losts, '-D', color='blue', label='Unique to Andromeda + Percolator at $10^{-2}$')
        ax.plot(fdr_ths, gains, '-X', color='green', label='Unique to '+tech+' + Percolator\n at given FDR (x-axis) ')
        plt.axhline(y=0, color='grey', linestyle='-', linewidth=0.75)
        plt.axhline(y=overlap_at_001, color='black', linestyle='--', linewidth=0.75)
        plt.axhline(y=1-overlap_at_001, color='grey', linestyle='--', linewidth=0.75)
        plt.axhline(y=1, color='grey', linestyle='-', linewidth=0.75)
        plt.axvline(x=0.01, color='grey', linestyle='--', linewidth=0.75)
        
        #Details
        ax.set_xscale('log', basex=10)
        plt.gca().invert_xaxis()
        plt.ylabel('Relative overlap with Andromeda at 1% FDR [in %]', fontsize=16)
        plt.xlabel('FDR threshold', fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        ax.legend(fontsize=16, frameon=True, loc='center left')
        plt.tight_layout()
        plt.savefig('gain_lost_plot.svg', format='svg')
        plt.show()
    
    return
