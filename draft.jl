# Lire les lignes du fichier
lines = readlines("scenarios/scenarios_h6.txt")

# Initialiser liste pour les valeurs non valides
invalid_values = []

# Parcourir chaque ligne
for (i, line) in enumerate(lines)
    nums = parse.(Int, split(line))  # convertir en entiers
    for (j, val) in enumerate(nums)
        if val != 0 && val != 1
            push!(invalid_values, (i, j, val))  # (ligne, colonne, valeur)
        end
    end
end

# Afficher les résultats
if !isempty(invalid_values)
    println("Valeurs différentes de 0 et 1 trouvées :")
    for (i, j, val) in invalid_values
        println("  À la ligne $i, colonne $j : $val")
    end
else
    println("Toutes les valeurs sont 0 ou 1.")
end