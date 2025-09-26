#
#
from numpy import abs, diff, linspace, interp
from pandas import DataFrame, concat, ExcelWriter, Series


def calc_L11(Relambda, ReLe):
    return 6.12556331 * (Relambda**-0.82218137) + 0.42350691
    # (Relambda**2 / ReLe) * (3 / 2) ** (1 / 2) / 20


def calc_U_max(U, A, D):
    """Calcula U max de acordo com a abertura mínima da válvula
    C é a abertura mínima
        w: Vazao kg/min
        A: Abertura de valvula 0 a 100
    """
    Ab = linspace(0.1, 1, 1000)
    y = linspace(0.0519 / 1000, D, 1000)
    Ab = linspace(0, 1, 1000)  # Abertura de 0 a 1
    y = linspace(0.0 / 1000, D, 1000)  # Altura de 0 a D
    h = interp(A / 100, Ab, y)  # altura interpolada para A
    U_max = U * D / h
    At_min = 0
    return U_max, At_min, h

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
    d43_exp = exp.get_D_caracteristico(dtg, sol["pbe_sol"].moc.xi_d, dff=[4, 3])
    diff_DTG = sum(abs(sol["pbe_sol"].N2Fv - dtg["freq_v"] / 100))
    # Função erro relativo para cada distribuição:
    if tipo == "psi":
        return sum((sol["pbe_sol"].N2Fv - dtg["freq_v"] / 100) ** 2) / sum(
            (dtg["freq_v"] / 100) ** 2
        )
    elif tipo == "Er":
        return (sol["pbe_sol"].moc.d43 - d43_exp) / (d43_exp)

    elif tipo == "D43":
        return (sol["pbe_sol"].moc.d43 - d43_exp) / d43_exp
    else:
        raise Exception(f"tipo not defined: {tipo}")


def get_properties(runs, solutions, IGNORAR_TESTES):

    psi = []
    Er = []
    run = []
    run_id = []
    nome= []
    N_emul = []
    N_esc = []
    marco = []
    method = []
    mv01 = []
    mv02 = []
    compare = []
    idict = []
    ch2o = []
    re = []
    we = []
    dt = []
    ts = []
    fdiss = []
    epsilon = []
    D43_r = []

    # xi = []
    # gamma = []

    j = 0
    for r in runs:
        sol = runs[r]
        print(f"Getting {r} simulation")
        for i in sol:  # Iterate over each entry
            if not isinstance(i, int):
                continue
            if sol[i]["N_escoam"] in IGNORAR_TESTES:
                continue
            if sol[i].get("opt_flag") is not None:
                print("Optmization case, better use another pos-process")
                break
            N_emul.append(int(sol[i]["exp"]["N_emul"]))
            N_esc.append(int(sol[i]["N_escoam"]))
            idict.append(i)
            run.append(r)
            nome.append(solutions[r])
            marco.append(sol[i]["marco"])
            if sol[i]["compares"] == [2, 3]:
                compare.append("E_ANM")
            if sol[i]["compares"] == [5, 6]:
                compare.append("E_Choke")
            mv01.append(sol[i]["exp"]["ANM"])
            mv02.append(sol[i]["exp"]["Choke"])
            dtg_a, d43_a, dtg_d, d43_d = get_dtg(sol["experiments"], sol[i])
            D43_r.append(d43_d / d43_a)
            we.append(dtg_a["We"].mean())
            re.append(sol[i]["exp"]["Re_MV01"])
            ch2o.append(sol[i]["exp"]["C_agua [%]"])
            epsilon.append(sol[i]["pbe_sol"].cp.epsilon)
            psi.append(calc_error(sol[i], sol["experiments"]))
            dt.append(diff(sol[i]["pbe_sol"].time).mean())
            ts.append(len(sol[i]["pbe_sol"].time))
            fdiss.append((sol[i]["pbe_sol"].fdiss))
            Er.append(calc_error(sol[i], sol["experiments"], "Er"))
            # xi.append(sol[i]["pbe_sol"].moc.xi)
            # gamma.append(sol[i]["pbe_sol"].moc.gamma)
            # beta = sol[i]["pbe_sol"].moc.beta
            # Q = sol[i]["pbe_sol"].moc.Q
            run_id.append(j)
        j += 1

    data = DataFrame(run, columns=["Parâmetro"])
    data = concat(
        [
            data,
            DataFrame(
                {
                    "Nome": nome,
                    "Erro": psi,
                    "Erro_Er": Er,
                    "run_id": run_id,
                    "N_emul": N_emul,
                    "N_esc": N_esc,
                    "Marco": marco,
                    "Compare": compare,
                    "MV01 [%]": mv01,
                    "MV02 [%]": mv02,
                    "epsilon": epsilon,
                    "Re": re,
                    "We": we,
                    "Red D43": D43_r,
                    "H2O [%]": ch2o,
                    "idict": idict,
                    "dt": dt,
                    "ts": ts,
                    "fdiss": fdiss,
                }
            ),
        ],
        axis=1,
    )
    data["Erro"] *= 100
    data["Erro_Er"] *= 100

    return data
