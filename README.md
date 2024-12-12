# Suppa
Controllare versione di Python
```bash
python --version
```
#Python 3.12.6

Creo ambiente python
```bash
python3 -m venv suppa_env
source suppa_env/bin/activate #attiva ambiente 
```

scarico code di [Suppav2.4](https://github.com/comprna/SUPPA/releases/tag/v2.4) e lo rendo eseguibile in ogni directory
```bash
alias suppa='python /Users/labbione/Desktop/SUPPA/SUPPA-2.4/suppa.py'
```
Verifico il funzionamento
```Python
suppa -h
```

## Generazione Eventi
```Python
suppa generateEvents -i /Volumes/Arianna/HeLa_MMC/Suppa2/gencode.v46.chr_patch_hapl_scaff.annotation.gtf -o /Volumes/Arianna/HeLa_MMC/Suppa2/New_suppa/Events/events.ioe -f ioe -e SE SS MX RI FL
```
mi crea un file per ogni evento, quindi copio l'intestazione dal primo file e aggiungo il contenuto di tutti i file, saltando l'intestazione
```Python
head -n 1 events.ioe_A3_strict.ioe > events.ioe
for file in events.ioe_*.ioe; do
    tail -n +2 "$file" >> events.ioe
done
```
Ho il file events.ioe
```
seqname	gene_id	event_id	alternative_transcripts	total_transcripts
chr1	ENSG00000237491.10	ENSG00000237491.10;A3:chr1:779092-803919:779092-803951:+	ENST00000655765.1,ENST00000670700.1,ENST00000657896.1,ENST00000658648.1,ENST00000656571.1,ENST00000669749.1,ENST00000666217.1,ENST00000665719.1	ENST00000656571.1,ENST00000666217.1,ENST00000657896.1,ENST00000655765.1,ENST00000412115.2,ENST00000669749.1,ENST00000658648.1,ENST00000665719.1,ENST00000670700.1
chr1	ENSG00000237491.10	ENSG00000237491.10;A3:chr1:779092-803919:779092-803922:+	ENST00000655765.1,ENST00000670700.1,ENST00000657896.1,ENST00000658648.1,ENST00000656571.1,ENST00000669749.1,ENST00000666217.1,ENST00000665719.1	ENST00000656571.1,ENST00000666217.1,ENST00000657896.1,ENST00000655765.1,ENST00000443772.2,ENST00000669749.1,ENST00000658648.1,ENST00000665719.1,ENST00000670700.1
chr1	ENSG00000237491.10	ENSG00000237491.10;A3:chr1:807323-809622:807323-809658:+	ENST00000412115.2,ENST00000670700.1,ENST00000657896.1,ENST00000658648.1	ENST00000656571.1,ENST00000657896.1,ENST00000412115.2,ENST00000669749.1,ENST00000658648.1,ENST00000588951.5,ENST00000665719.1,ENST00000670700.1
[...]
```
> [!IMPORTANT]
> Controllare che tutti i file tpm abbiano valori di espressione numerici e tutti con le stesse colonne, creare lo script check_tpm_values.py
```Python
import pandas as pd
import glob

def check_file_format(filename):
    try:
        # Legge il file
        df = pd.read_csv(filename, sep='\t', header=None)
        
        # Controlla il numero di colonne
        num_columns = len(df.columns)
        if num_columns != 2:
            print(f'ATTENZIONE: Il file {filename} ha {num_columns} colonne invece delle 2 attese')
            print(f'Prime righe del file {filename}:')
            print(df.head())
            return None
            
        # Se il formato è corretto, procedi con l'analisi
        df.columns = ['Transcript_ID', 'Expression_Value']
        
        # Verifica i valori numerici
        non_numeric = pd.to_numeric(df['Expression_Value'], errors='coerce').isna()
        
        if non_numeric.any():
            print(f'Nel file {filename} ci sono valori non numerici:')
            print(df[non_numeric])
        else:
            print(f'Nel file {filename} tutti i valori di espressione sono numerici.')
            
        return df
        
    except Exception as e:
        print(f'Errore nel processare il file {filename}:')
        print(str(e))
        return None

# Processa tutti i file .tpm
for filename in sorted(glob.glob('*.tpm')):
    print(f'\nAnalizzando {filename}...')
    df = check_file_format(filename)
    if df is not None:
        print(f'Numero di righe: {len(df)}')
        print(f'Range dei valori: {df["Expression_Value"].min()} - {df["Expression_Value"].max()}')
```
```Python
python check_tpm_values.py
```




suppa psiPerEvent -i /Volumes/Arianna/HeLa_MMC/Suppa2/New_suppa/filtered/filtered_events.ioe -e /Volumes/Arianna/HeLa_MMC/Suppa2/tpm_files/A_all.tpm -o /Volumes/Arianna/HeLa_MMC/Suppa2/New_suppa/Psi/A_psi
