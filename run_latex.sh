#!/bin/bash

BUILD_DIR="build"
LOG_FILE="dump.txt"
mkdir -p "$BUILD_DIR"
HAS_ERRORS=0

# Iteriere über alle .tex Dateien im aktuellen Ordner
for tex_file in *.tex; do
    # Prüfen, ob überhaupt .tex Dateien gefunden wurden
    [ -e "$tex_file" ] || continue

    # --- AUSNAHMEN ---
    if [[ "$tex_file" == tensor_definition_order*.tex ]]; then
        continue
    fi
    # --------------------
    base_name=$(basename "$tex_file" .tex)
    echo "Kompiliere $tex_file..."
    # Kompilieren: Ausgabe umleiten (stdout und stderr) in dump.txt
    # -output-directory sorgt dafür, dass Hilfsdateien nicht im Hauptordner liegen
    # -interaction=nonstopmode verhindert, dass LaTeX bei Fehlern anhält
    if pdflatex -interaction=nonstopmode -output-directory="$BUILD_DIR" "$tex_file" > "$LOG_FILE" 2>&1; then
        echo "  -> Successfully compiled."
    else
        echo "  -> ERROR for $tex_file (see $LOG_FILE)."
        HAS_ERRORS=1
    fi
done

# Nach dem Durchlauf prüfen, ob alles erfolgreich war
if [ $HAS_ERRORS -eq 0 ]; then
    echo "All Files compiled. I delete the log..."
    rm "$LOG_FILE"
else
    echo "Some File Compilations failed. I keep '$LOG_FILE'."
fi
echo "Done."