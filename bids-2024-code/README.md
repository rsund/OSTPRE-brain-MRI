# Tiedostot
- tags.csv 
    - Kaikki metadatat dicomeista

- vanhat.csv
    - Metadata dicomeista, jotka käännettiin bidsiin. Vanha bids-polku löytyy sarakkeesta ['bidsPath']
    - Sisältää cat:sta saadun laatuarvosanan ja prosentin ["grade","percentage"]

- uudet.csv
    - Metadata dicomeista, jotka saatiin lisäämällä tarkistettavia sarakkeita ['SeriesDescription','PulseSequenceName'] alunperin tarkastettiin vain ['ProtocolName']

- uudet-uniikit.csv
    - vanhat dicomit poistettu

- "new-bids"-kansio
    - sisältää yhden T1w-kuvan per sessio 
    - jos useampia T1w:ta per sessio, niin ne ylikirjoittuvat (pitäisikö kaikki tallentaa?)
        - esimerkiksi näistä kuvista tallentui vain yksi:
            - sub-00532/ses-01/ses-1_t1_mprage_sag_p2_iso_1.0_20180313122235_12
            - sub-00532/ses-01/ses-1_t1_mprage_sag_p2_iso_1.0_20180313122235_8
            - sub-00532/ses-01/ses-1_t1_mprage_sag_p2_iso_1.0_MPR_Tra_20180313122235_14


# Skriptit
- filter-metadata.py
    - luo uudet.csv:n

- find-unique.py
    - hakee kuvat, jotka eivät ole olleet analyysissä mukana
    - luo uudet-uniikit.csv:n

- create-bids-folders.py
    - luo bids-tiedostot niftien perusteella
        - 2663 kuvan konvertointi tuottaa 1201 bids-tiedostoa
            - kuvia katoaa mm. ylikirjoittaisen vuoksi, puuttuvien tietojen vuoksi (ikä / sukupuoli) ja muista syistä, joita en ole tutkinut tarkemmin

- plot-size.py
    - histogrammi tiedostokoon ja kuvamäärän perusteella, vertailee uusia ja vanhoja kuvia




