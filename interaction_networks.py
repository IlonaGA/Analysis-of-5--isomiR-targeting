import pandas as pd
import numpy as np
import os
import py4cytoscape as p4c
p4c.cytoscape_ping()
p4c.cytoscape_version_info()

projects = [file.split('_')[0] for file in os.listdir('directory')]  #Directory with raw and summary files (obtained by 1_create_raw_and_summary.py).
projects = set(projects)

databases = ['miRDB_and_TargetScan']

for database in databases:  
    projects = [file.split('_')[0] for file in os.listdir('TMM-RPM')]
    projects = set(projects)
    #projects = ['TCGA-COAD']
    for project in projects:
        try:
            data = pd.read_csv(f'TMM-RPM/{project}_tumor_raw.tsv', sep='\t', index_col=0)
            highly_expressed = pd.read_csv(f'highly_expressed_isomiRs/{project}_tumor.tsv', sep='\t', index_col=0) #List of highly expressed miRNAs.


            data = data.loc[set(data.index) & set(highly_expressed.index)]
            data.reset_index(inplace=True)
            data = data[data['corr'] < -0.3]

            isomirs = list(set(data['isomiR']))
            genes = list(set(data['gene']))
            nodes = pd.DataFrame(data={'id': isomirs + genes, 'type': ['isomiR'] * len(isomirs) + ['gene'] * len(genes)})
            edges = pd.DataFrame(data={'source': data['isomiR'], 'target': data['gene'], 'weight': [1] * len(data)})
            node_degree = pd.concat([edges['source'].value_counts(), edges['target'].value_counts()]) 
            node_degree = np.log(node_degree * 10) * 20
            node_degree.rename('degree', inplace=True)
            nodes = nodes.set_index('id').join(node_degree).reset_index()
            p4c.create_network_from_data_frames(nodes, edges, title=project, collection="DataFrame Example")

            p4c.set_visual_style('Nested Network Style')
            #p4c.set_node_color_mapping(**gen_node_color_map('type', mapping_type='d', style_name=style_nameb)) 
            p4c.set_node_color_mapping('type', table_column_values=['gene', 'isomiR'], colors=['#b6d8f2', '#7DAEF3'], mapping_type='d', style_name='Nested Network Style')
            #p4c.set_node_size_mapping('id', table_column_values=nodes['degree'].unique().tolist(), sizes=[60], mapping_type='d')#, style_name=style_name)
            p4c.set_node_size_mapping('id', nodes['id'].tolist(), nodes['degree'].tolist(), mapping_type='d', style_name='Nested Network Style')
            p4c.layout_network('force-directed-cl')
            p4c.set_node_shape_default('ELLIPSE')
            p4c.export_image(f'{project}', type='PDF', overwrite_file=True)
            p4c.sandbox.sandbox_get_from(f'{project}.PDF', f'/output/directory/{project}.pdf')
        except:
            print(project)
    



