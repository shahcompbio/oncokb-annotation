---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.16.1
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Settings

```python
import os
import pandas as pd
```

# Data

```python
patient = '107C1'
vartype = 'indel'
root_dir = '/data1/shahs3/users/chois7/projects/colibactin/signatures'
tables_dir = f'{root_dir}/{patient}/tables'
```

## Make annovar inputs

```python
patients = ['107C1', '136PC2']
vartypes = ['snv', 'indel']
comp_sets = [
    'ClbPpos - dClbP - Ctrl - W0',
    'dClbP - ClbPpos - Ctrl - W0',
    'Ctrl - ClbPpos - dClbP - W0',
    'W0 - ClbPpos - dClbP - Ctrl',
]
for patient in patients:
    for vartype in vartypes:
        for comp_set in comp_sets:
            var_path = f'{tables_dir}/unique_dfs.{comp_set}.{vartype}.all.tsv'
            var = pd.read_table(var_path)
            with open(f'../results/avinput/{patient}.{comp_set.replace(" ", "")}.{vartype}.avinput', 'w') as out:
                for rix, row in var.iterrows():
                    chrom = row['chrom']
                    pos = row['pos']
                    ref = row['ref']
                    alt = row['alt']
            
                    if vartype == 'snv':
                        start = end = pos
                    elif vartype == 'indel':
                        is_ins = len(alt) > len(ref)
                        if is_ins:
                            assert len(ref) == 1
                            start = end = pos
                            ref = '-'
                            alt = alt[1:]
                        else:
                            assert len(alt) == 1 and len(ref) > 1
                            start = pos+1
                            end = pos + len(ref)
                            ref = ref[1:]
                            alt = '-'
                    else:
                        raise ValueError(vartype)
                        
                    field = [chrom, start, end, ref, alt] # 1       948921  948921  T       C
                    line = '\t'.join([str(s) for s in field]) + '\n'
                    out.write(line)
```

# Oncokb

```python
import requests
```

```python
def get_oncokb_response(gene, protein_alteration, api_key):
    genome = 'GRCh38'

    headers = {
        'accept': 'application/json',
        'Authorization': f'Bearer {api_key}',
    }
    
    params = {
        'hugoSymbol': gene,
        'alteration': protein_alteration,
        'referenceGenome': genome,
    }
    
    api_url = 'https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange'
    response = requests.get(api_url, params=params, headers=headers)
    
    if response.ok:
        data = response.json()
        return data
```

```python
def test_get_oncokb_response():
    api_key = str(os.environ['ONCOKB_API_KEY']).rstrip()
    gene = 'BRAF'
    protein_alteration = 'p.V600E'
    data = get_oncokb_response(gene, protein_alteration, api_key)
    assert data['oncogenic'] == 'Oncogenic', data

test_get_oncokb_response()
```

```python
import os
```

```python
os.environ['ONCOKB_API_KEY'] = open('/home/chois7/keys/.oncokb_apikey').read().rstrip()
```

```python
api_key = str(os.environ['ONCOKB_API_KEY']).rstrip()
for patient in patients:
    for comp_set in comp_sets:
        for vartype in vartypes:
            df = pd.DataFrame(columns=['gene', 'refseq', 'alteration', 'oncogenic'])
            in_path = f'../results/avinput/{patient}.{comp_set.replace(" ", "")}.{vartype}.avinput.exonic_variant_function'
            out_path = f'../results/oncokb/{patient}.{comp_set.replace(" ", "")}.{vartype}.oncokb.tsv'
            for line in open(in_path, 'r'):
                field = line.strip().split('\t')
                annots = field[2].strip(',').split(',')
                for annot in annots:
                    gene, gene_id, exon, cdna_alt, protein_alt = annot.split(':')
                    oncokb = get_oncokb_response(gene, protein_alt, api_key)
                    row = [gene, gene_id, protein_alt, oncokb['oncogenic']]
                    df.loc[df.shape[0]] = row
            df.to_csv(out_path, sep='\t', index=False)
```

```python
df['oncogenic'].value_counts()
```

```python
gene, protein_alt
```

```python
line.strip().split('\t')
```
