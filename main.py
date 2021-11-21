import numpy as np
import streamlit as st

st.write("Cálculo de los parámetros GSI, D y mi del criterio de falla de Hoek-Brown a partir de la ecuación potencial normalizada de la envolvente de falla en el plano τ vs. σ.")

A = st.number_input('Valor de "A" normalizado:',value=0.16113, format="%.10f")
B = st.number_input('Valor de "B" normalizado:',value=0.71346, format="%.10f")
stn = st.number_input('Valor del esfuerzo de tensión normalizado "𝜎𝑡n"',value=-2.53E-04, format="%.10f")
sci = st.number_input('Valor de resistencia a la compresión inconfinada de la roca intacta "𝜎ci [kPa]"',value=99600.0, format="%.10f")


def evaluate(A, B, sci, stn):
    
    n_dots = 12500
    n = np.linspace(0.1, 60, n_dots)
    sn = n / sci

    tn = A * (sn - stn) ** B
    t = tn * sci
    
    beta = np.arctan(A * B * (sn - stn) ** (B - 1))
    sprima = ((2 * np.sin(beta) * (1 + np.sin(beta))) / (np.cos(beta)) ** 2) + 1
    s1 = n + t * sprima ** 0.5
    s3 = n - t / sprima ** 0.5

    GSI = np.linspace(1, 100, n_dots)
    
    result = []
    for i, gsi in enumerate(GSI):        
        a = 0.5 + 1 / 6 * (np.exp(-gsi / 15) - np.exp(-20 / 3))
        mb = (sprima[0] - 1) / (a * ((s1[0] - s3[0]) / sci ) ** ((a - 1) / a))
        s = ((s1[0] - s3[0]) / sci ) ** (1 / a) - mb * s3[0] / sci 
        if 0 < s:
            D = ((100 - gsi) / (3 * np.log(s))) + 3
            if 0 < D and D <= 1:
                mi = mb / np.exp((gsi - 100) / (28 - 14 * D))
                s1_teo = s3 + sci * (mb * s3 / sci + s) ** a
                if 0 < mi and mi <= 100:
                    result.append(
                        {
                            "D": D,
                            "mi": mi,
                            "GSI": gsi,
                            "s1": s1_teo.tolist(),
                            "res": np.mean(abs(s1_teo - s1)),
                            }
                        )

    min_res = min(result, key=lambda x: x["res"])
    return min_res


result = evaluate(A, B, sci, stn)

st.title("Resultados")
st.write(f'GSI: {result["GSI"]:.0f}')
st.write(f'D: {result["D"]:.2f}')
st.write(f'mi: {result["mi"]:.2f}')
st.write(f'error: {result["res"]:.4f}')

st.title("Referencias")
st.write('Hoek E., Brown E.T. (1980) “Empirical strength criterion for rock masses”. J. Geotech. Engng Div., ASCE 106(GT9), 1013-1035')
st.write('Kumar, P. 1999. Shear Failure Envelope of Hoek-Brown Criterion for Rockmass')
       
