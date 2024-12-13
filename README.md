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

scarico code di [Suppav2.4](https://github.com/comprna/SUPPA/releases/tag/v2.4) e lo rendo eseguibile in ogni directory (ogni volta che chiudo o riapro il terminale)
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
## Files TPM
> [!IMPORTANT]
> Controllare che tutti i file tpm abbiano valori di espressione numerici e tutti con le stesse colonne, creare lo script `check_tpm_values.py`
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
Genero file TPM unito per ogni replica
```Python
#Replica A
suppa joinFiles -f tpm -i A1.tpm A2.tpm A3.tpm A4.tpm -o A_all
#Replica B
suppa joinFiles -f tpm -i B1.tpm B2.tpm B3.tpm B4.tpm -o B_all
#Replica c
suppa joinFiles -f tpm -i C1.tpm C2.tpm C3.tpm C4.tpm -o C_all
```


Dobbiamo verificare che i trascritti espressi siano consistenti tra tutti i campioni (A, B, e C) e creare un file `filtered_events.ioe` che contenga solo gli eventi con trascritti presenti in tutti i file di espressione.\
Creo uno script `analyze_all_replicates.py`
```Python
import pandas as pd
import os
from datetime import datetime

def analyze_transcripts_all_replicates():
    print(f"Analisi iniziata il: {datetime.utcnow()}")
    
    # Percorsi dei file
    events_file = '/Volumes/Arianna/HeLa_MMC/Suppa2/New_suppa/Events/events.ioe'
    expression_files = {
        'A': '/Volumes/Arianna/HeLa_MMC/Suppa2/tpm_files/A_all.tpm',
        'B': '/Volumes/Arianna/HeLa_MMC/Suppa2/tpm_files/B_all.tpm',
        'C': '/Volumes/Arianna/HeLa_MMC/Suppa2/tpm_files/C_all.tpm'
    }
    output_dir = '/Volumes/Arianna/HeLa_MMC/Suppa2/New_suppa/filtered'
    
    # Crea la directory di output se non esiste
    os.makedirs(output_dir, exist_ok=True)
    
    print("\nVerificando i file di input...")
    for name, path in expression_files.items():
        if os.path.exists(path):
            print(f"File {name}_all.tpm trovato")
        else:
            print(f"ERRORE: File {name}_all.tpm non trovato in {path}")
            return
            
    if os.path.exists(events_file):
        print("File events.ioe trovato")
    else:
        print(f"ERRORE: File events.ioe non trovato in {events_file}")
        return

    print("\nLeggendo i file...")
    
    # Leggi il file events.ioe
    events = pd.read_csv(events_file, sep='\t')
    print(f"Eventi totali nel file originale: {len(events)}")
    
    # Raccogli tutti i trascritti dal file events.ioe
    event_transcripts = set()
    for idx, row in events.iterrows():
        event_transcripts.update(row['alternative_transcripts'].split(','))
        event_transcripts.update(row['total_transcripts'].split(','))
    
    print(f"Trascritti totali negli eventi: {len(event_transcripts)}")
    
    # Leggi i file di espressione e trova i trascritti comuni
    expression_transcripts = {}
    for replicate, file_path in expression_files.items():
        expression = pd.read_csv(file_path, sep='\t', index_col=0)
        expression_transcripts[replicate] = set(expression.index)
        
    # Trova i trascritti presenti in tutti i file di espressione
    common_transcripts = set.intersection(*expression_transcripts.values())
    
    # Statistiche dettagliate
    print("\nSTATISTICHE DETTAGLIATE PER REPLICATO:")
    for replicate, transcripts in expression_transcripts.items():
        print(f"\nReplicato {replicate}:")
        print(f"Numero totale di trascritti: {len(transcripts)}")
        missing = event_transcripts - transcripts
        unique = transcripts - set.union(*(expression_transcripts[r] for r in expression_transcripts if r != replicate))
        print(f"Trascritti mancanti rispetto agli eventi: {len(missing)}")
        print(f"Trascritti unici in questo replicato: {len(unique)}")
    
    print(f"\nTrascritti comuni a tutti i replicati: {len(common_transcripts)}")
    print(f"Percentuale di trascritti comuni: {len(common_transcripts)*100/len(event_transcripts):.2f}%")
    
    # Filtra gli eventi
    print("\nFiltraggio eventi...")
    def check_transcripts(row, available_transcripts):
        alt_trans = set(row['alternative_transcripts'].split(','))
        total_trans = set(row['total_transcripts'].split(','))
        return all(t in available_transcripts for t in alt_trans.union(total_trans))
    
    events['valid'] = events.apply(lambda row: check_transcripts(row, common_transcripts), axis=1)
    
    # Salva gli eventi validi
    valid_events = events[events['valid']].drop(columns=['valid'])
    output_file = os.path.join(output_dir, 'filtered_events_all_replicates.ioe')
    valid_events.to_csv(output_file, sep='\t', index=False)
    
    # Analisi per tipo di evento
    print("\nANALISI PER TIPO DI EVENTO:")
    events['event_type'] = [e.split(';')[1].split(':')[0] for e in events['event_id']]
    valid_events['event_type'] = [e.split(';')[1].split(':')[0] for e in valid_events['event_id']]
    
    print("\nEventi originali per tipo:")
    print(events['event_type'].value_counts())
    print("\nEventi filtrati per tipo:")
    print(valid_events['event_type'].value_counts())
    
    print(f"\nAnalisi completata il: {datetime.utcnow()}")
    
    return output_file

if __name__ == "__main__":
    filtered_file = analyze_transcripts_all_replicates()
    if filtered_file:
        print(f"\nFile degli eventi filtrato creato: {filtered_file}")
        print("\nComandi per eseguire SUPPA:")
        for replicate in ['A', 'B', 'C']:
            print(f"\nPer il replicato {replicate}:")
            print(f"suppa psiPerEvent -i {filtered_file} \\\n-e /Volumes/Arianna/HeLa_MMC/Suppa2/tpm_files/{replicate}_all.tpm \\\n-o /Volumes/Arianna/HeLa_MMC/Suppa2/New_suppa/Psi/{replicate}_psi")
```
```Python
python analyze_all_replicates.py
```
Risultato:\
Eventi totali nel file originale: 327628\
Trascritti totali negli eventi: 199209\

### STATISTICHE DETTAGLIATE PER REPLICATO:
\
Replicato A:\
Numero totale di trascritti: 253181\
Trascritti mancanti rispetto agli eventi: 16121\
Trascritti unici in questo replicato: 0\
\
Replicato B:\
Numero totale di trascritti: 253181\
Trascritti mancanti rispetto agli eventi: 16121\
Trascritti unici in questo replicato: 0\
\
Replicato C:\
Numero totale di trascritti: 253181\
Trascritti mancanti rispetto agli eventi: 16121\
Trascritti unici in questo replicato: 0\
\
Trascritti comuni a tutti i replicati: 253181\
Percentuale di trascritti comuni: 127.09%\

### ANALISI PER TIPO DI EVENTO:

Eventi originali per tipo:\
event_type\
AF    139285\
SE     61582\
AL     52803\
A3     27004\
A5     23789\
RI     12011\
MX     11154\
Name: count, dtype: int64\
\
Eventi filtrati per tipo:\
event_type\
AF    130713\
SE     56795\
AL     45498\
A3     24248\
A5     21695\
RI     10548\
MX     10434

- Tutti e tre i replicati (A, B, C) hanno esattamente lo stesso numero di trascritti: 253,181
- Hanno tutti lo stesso numero di trascritti mancanti: 16,121
- Non ci sono trascritti unici in nessun replicato (tutti 0)
- I trascritti sono perfettamente consistenti tra i replicati

Statistiche degli eventi: \
Eventi originali totali: 327,628 \
Eventi dopo il filtraggio: \
AF: da 139,285 a 130,713 (-8,572) \
SE: da 61,582 a 56,795 (-4,787) \
AL: da 52,803 a 45,498 (-7,305) \
A3: da 27,004 a 24,248 (-2,756)\
A5: da 23,789 a 21,695 (-2,094) \
RI: da 12,011 a 10,548 (-1,463) \
MX: da 11,154 a 10,434 (-720) \
Quindi la qualità dei dati è buona è il processo di quantificazione è stato eseguito in modo coerente.
