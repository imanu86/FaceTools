# Checklist test manuale GUI (da eseguire con il build v11 installato)
Usa un progetto di prova (il frame sintetico va bene). Annota bug nella colonna Bug.

| # | Passo | Atteso | Bug |
|---|---|---|---|
| 1 | **Cut and extend – U/Y**: seleziona guide su un ramo | Si muove solo la regione collegata al ramo; il resto fermo, nessun fallback | ☐ |
| 2 | **T**: due bracci della stessa barra superiore | Si muove l'intera barra collegata | ☐ |
| 3 | Guide ambigue / fisse / non associate | Rifiuto con messaggio, nessuna scelta silenziosa | ☐ |
| 4 | Cambia sorgente/posa/piano/maschere/ambito, nuovo piano, clear | Il consenso e il lavoro vecchi si invalidano | ☐ |
| 5 | Solo cambio distanza nella stessa regione | Consenso mantenuto; 0 = no-op; negativo rifiutato | ☐ |
| 6 | Prima visualizzazione prima della conferma | Tutta la regione arancione è disegnata e ispezionabile | ☐ |
| 7 | Apply, poi Undo e Redo | Posa sorgente, resto fermo e altre parti intatti; lavoro rifiutato non cambia nulla | ☐ |
| 8 | Esporta (strict export): annulla / errore / ok | Modale non stale, file di output corretto | ☐ |
| 9 | Salva progetto 3MF, riapri, slice | Parti identiche, slicing ok | ☐ |
| 10 | Estrusione faccia: click / Shift-click | Selezione facce stesso part; normali opposte rifiutate | ☐ |
| 11 | Trascina distanza | Anteprima e badge si aggiornano durante il drag | ☐ |
| 12 | Incrementi 1 mm relativi; input numerico | Valore frazionario iniziale preservato | ☐ |
| 13 | Orbit/pan con facce selezionate; Esc/focus perso | Selezione mantenuta; pulizia sicura | ☐ |
| 14 | Nuova parte / Union / Subtract + Undo/Redo singolo | Parti non correlate e fori intatti | ☐ |
| 15 | Chiudi/ricarica/undo durante calcolo | Nessun risultato stale applicato | ☐ |
| 16 | Parte negativa singola: attiva Faces, click | Cutter isolato; negative miste/multiple rifiutate | ☐ |
| 17 | Negativa: NewPart negativo, Union/Subtract, salva/riapri, slice | Cavità corretta, parti positive intatte | ☐ |
| 18 | Posizionamento per coincidenza di punti | Parte allineata; Undo ripristina | ☐ |
| 19 | Smart fill by angle | Selezione coerente con l'angolo | ☐ |

Tutte "solo verifica" (build già installato). Se un test fallisce → correzione codice + build C++ (ore).
