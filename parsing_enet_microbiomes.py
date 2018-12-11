import pandas as pd


ref_data= pd.read_csv('taxonomy_7_levels.txt',header=None,sep='\t')
ref_data.columns = ['taxid','names']
data = pd.read_csv('enet_microbiomes.txt',header=None)
data.columns=['taxid']


out = pd.merge(data,ref_data,on='taxid')
bname_output = pd.DataFrame()
for index,row in data.iterrows():
    string=row['taxid']
    if '+' not in string:
        array=string.split(' ')
        output=ref_data.loc[ref_data['taxid'].isin(array)]
        output.insert(0,'index_id',[str(index)]*len(output.index))
        bname_output=bname_output.append(output)

    else:
        array = string.split('+')
        output = ref_data.loc[ref_data['taxid'].isin(array)]
        output.insert(0, 'index_id', [str(index)] * len(output.index))
        bname_output = bname_output.append(output)


divisons=bname_output['names'].str.split(';',n=6,expand=True)
bname_output["Domain"]=divisons[0]
bname_output["Kingdom"]=divisons[1]
bname_output["Phylum"]=divisons[2]
bname_output["Class"]=divisons[3]
bname_output["Family"]=divisons[4]
bname_output["Genus"]=divisons[5]
bname_output["Species"]=divisons[6]
bname_output.drop(columns=['names'],inplace=True)
bname_output.groupby(['index_id','Domain'])
lvl0 = bname_output.index_id.values
lvl1 = bname_output.Domain.values
lvl2 = bname_output.Kingdom.values
lvl3 = bname_output.Phylum.values
lvl4 = bname_output.Class.values
lvl5 = bname_output.Family.values
lvl6 = bname_output.Genus.values
lvl7 = bname_output.Species.values
lvl8 = bname_output.taxid.values
midx=pd.MultiIndex.from_arrays([lvl0,lvl1,lvl2,lvl3,lvl4,lvl5,lvl6,lvl7,lvl8],
                               names=['index_id','Domain','Kingodm','Phylum','Class','Family','Genus','Species','taxid'])
col=['Filler']
df = pd.DataFrame('-',midx,col)




df.to_excel('elastic_net_microbiomes.xlsx')







