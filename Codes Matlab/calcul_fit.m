% Fonction de calcul du fit
function fit_value = calcul_fit(mesures, predictions)
    fit_value = 100 * (1 - norm(mesures - predictions) / norm(mesures - mean(mesures)));
end