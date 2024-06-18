#
#
from numpy import abs


#
#
def get_dtg(exp, sol):
    dtg_a = exp.get_DTG(
        teste=sol["N_escoam"], marco=sol["marco"], ID=sol["compares"][0]
    )
    dtg_d = exp.get_DTG(
        teste=sol["N_escoam"], marco=sol["marco"], ID=sol["compares"][1]
    )
    d43_a = exp.get_D_caracteristico(dtg_a, sol["pbe_sol"].moc.xi_d, dff=[4, 3])
    d43_d = exp.get_D_caracteristico(dtg_d, sol["pbe_sol"].moc.xi_d, dff=[4, 3])
    return dtg_a, d43_a, dtg_d, d43_d


def calc_error(sol, exp, tipo="psi"):
    """
    Args:
        sol (_type_):
        exp (_type_):
        tipo (str): Er or psi. Defaults to "psi".
    """
    dtg = exp.get_DTG(teste=sol["N_escoam"], marco=sol["marco"], ID=sol["compares"][1])
    d43 = exp.get_D_caracteristico(dtg, sol["pbe_sol"].moc.xi_d, dff=[4, 3])
    diff_DTG = sum(abs(sol["pbe_sol"].N2Fv - dtg["freq_v"] / 100))
    # Função erro relativo para cada distribuição:
    erro_D43 = abs(sol["pbe_sol"].moc.d43 - d43) / d43
    if tipo == "psi":
        return sum((sol["pbe_sol"].N2Fv - dtg["freq_v"] / 100) ** 2) / sum(
            (dtg["freq_v"] / 100) ** 2
        )
    elif tipo == "Er":
        return sum(sol["pbe_sol"].N2Fv - dtg["freq_v"] / 100) / sum(dtg["freq_v"] / 100)

    elif tipo == "D43":
        pass
    else:
        raise Exception(f"tipo not defined: {tipo}")
