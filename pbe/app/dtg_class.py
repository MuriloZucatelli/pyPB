import pandas as pd

# from .classes import FLOW
from numpy import arange, pi, array, interp, polyfit, polyval
from pbe.setup.system import Domain, DispersePhase, ContinuosPhase, FLOW, Prop, Fluid2
from pbe.solvers.moc import MOCSolution
from pbe.models import breakup, coalescence
from os.path import join
from pathlib import Path
from os import environ
from typing import Union, List
from sys import path
from Circuito.setup import data_def as dd_c
from Bettersizer.setup import data_def as dd_b

# TODO: implementar uma classe que ja traz todos os parâmetros constantes
#       do caso, pra não precisar, por exemplo, de obter o objeto
#       continuousphase pra depois calcular epsilon. Além de não precisar
#       salvar em self.D ou L parâmetros que não se encaixam nem na classe solução
#       nem em domain ou phase continua ou dispersa


class DTGSolution:
    def __init__(
        self,
        M=10,  # number of classes
        time=None,
        exp: pd.DateFrame = None,
        data=None,
        model_parameters=None,
    ):
        self.D = 0.02095  # [m] pipe diameter diameter
        self.L = 2.5  # [m] length diameter

        # oil
        self.cp = ContinuosPhase(
            name=data.tableprop.oil.name,
            rho=exp[data.tableprop.oil.name + " [kg/m³]"],  # [kg/m³]
            mu=exp[data.tableprop.oil.name + " [Pa.s]"],
        )  # [P = kg * m^-1 s^-1]  dynamic viscosity

        self.cp.epsilon = calc_epsilon(U, self.D, self.cp.mu, self.cp.rho)

        # water solution
        self.dp = DispersePhase(
            name=data.tableprop.water.name,
            phi=exp["C_agua [%]"],
            rho=exp[data.tableprop.water.name + " [kg/m³]"],  # [kg/m3]
            mu=exp[data.tableprop.water.name + " [Pa.s]"],
            sigma=exp["sigma [N/m]"],  # [P = kg * m^-1 s^-1]
        )

        self.domain = Domain(theta=theta, V=pi * self.L * (self.D / 2) ** 2, M=M)

        self.set_parameters(model_parameters)
        self.nf0()
        vmin = None

        # Função distribuição de gotas filhas
        beta = breakup.DDSD.coulaloglou_beta

        # Função frequencia de quebra
        g = breakup.breakupModels(
            name="coulaloglou", C=self.C, domain=self.domain, cp=self.cp, dp=self.dp
        ).gamma
        # Função frequencia de Coalescencia
        Qf = coalescence.coalescenceModels(
            name="coulaloglou", C=self.C, domain=self.domain, cp=self.cp, dp=self.dp
        ).Q

        # Initial probability density function distribuition
        # A0 = DSD.analitico(dp=self.dp).A0

        self.moc = MOCSolution(
            M,
            time,
            self.dp.v_max / M,
            n0=self.n0,
            xi0=vmin,
            beta=beta,
            gamma=g,
            Q=Qf,
            theta=self.domain.theta,
        )

    def set_parameters(self, model_parameters):
        if model_parameters is None:
            self.C = [0.4, 0.08, 2.8, 1.83e13]
        else:
            self.C = model_parameters

    def n0(self, v):
        return 0 * v

    def nf0(self):
        self.nf0 = self.domain.V / self.domain.theta

    @property
    def pbe_phi(self):
        return self.moc.total_volume / self.domain.V


# Function used by class
def calc_epsilon(U: float, D: float, mu: float, rho: float):
    """Calcula as propriedades turbulentas

    Args:
        U (float): average flow velocity
        D (float):
        mu (float):
        rho (float):

    Returns:
        float: epsilon
    """
    Re = U * D / mu * rho
    TI = 0.16 * Re ** (-1.0 / 8.0)
    u_rms = U * TI
    k = 3.0 / 2.0 * u_rms**2
    L_t = 0.038 * D
    epsilon = 0.09 * k ** (3.0 / 2.0) / L_t

    return epsilon


"""
    Classes to import flow sensor and DSD measures
"""


class Import_flow_DSD2:
    def __init__(
        self,
        folder: str,
        teste: dict = None,
    ) -> None:
        """Importa os dados da planta e a DTG

        Args:
            folder (str): Pasta com os dados
            teste (List[int], optional): Teste. Defaults to None.
        """
        self.folder = folder
        self.teste = list(teste.keys())

        if isinstance(teste, int):
            self.teste = [teste]
        self.marco = teste  # Rastreia teste e marco dos testes

        f1 = join(folder, "flow.xlsx")
        # Le os dados externos (Tipo,Weber,Concentração,...)
        df1 = pd.read_excel(f1, sheet_name="flow_ext")
        # Le os dados de media e desvio padrao (Tipo,Weber,Concentração,...)
        df2 = pd.read_excel(f1, sheet_name="flow_media")
        df3 = pd.read_excel(f1, sheet_name="flow_std")
        # le os dados do circuito (P,T,V,...)
        df4 = pd.read_excel(f1, sheet_name="flow")

        self.flow = FLOW(flow=df4, flow_ext=df1, flow_mean=df2, flow_std=df3)

        f2 = join(folder, "DTG.xlsx")
        self.DTG = pd.read_excel(f2, sheet_name="DTG")

        if teste is not None:
            self.get_teste()
            self.get_marco()

    def get_teste(self):
        """Extrai o/os testes para analise de BP"""

        def grupa_t(df: pd.DataFrame, key: str):
            dfn = pd.DataFrame()
            dfg = df.groupby(key)
            testes = list(set(dfg.groups) & set(self.teste))
            for t in testes:
                dfn = pd.concat([dfn, dfg.get_group(t)], ignore_index=True)
            return dfn

        x = self.flow
        col_b = dd_b.DTG_df_name["numero_escoamento"]  # = "N_escoam"
        col_c = dd_c.flow_ext_df_name["numero_escoamento"]  # = "N_escoam"
        try:
            self.flow.flow = grupa_t(x.flow, col_c)
            self.flow.flow_ext = grupa_t(x.flow_ext, col_c)
            self.flow.flow_mean = grupa_t(x.flow_mean, col_c)
            self.flow.flow_std = grupa_t(x.flow_std, col_c)
            self.DTG = grupa_t(self.DTG, col_b)
        except Exception:
            pass
        return self

    def get_marco(self):
        """Extrai o/os marcos para analise de BP"""

        def grupa_m(df: pd.DataFrame, key: str):
            dfn = pd.DataFrame()
            N_esc = dd_c.flow_ext_df_name["numero_escoamento"]
            for t in self.teste:
                dfg = df.groupby(by=N_esc).get_group(t)
                dfgm = dfg.groupby(key)
                marcos = list(set(dfgm.groups) & set(self.marco[t]))
                for m in marcos:
                    dfn = pd.concat([dfn, dfgm.get_group(m)], ignore_index=True)
            return dfn

        def ajuste(df: pd.DataFrame, fe: pd.DataFrame):
            dfn = pd.DataFrame()
            ti = dd_c.flow_ext_df_name["tempo_inicial"]
            tf = dd_c.flow_ext_df_name["tempo_final"]
            faixa = pd.concat([fe[ti], fe[tf]], axis=1)
            for i in range(len(faixa)):
                idi = (
                    (df["Timestamp"] - faixa.iloc[i].to_numpy().min())
                    .abs()
                    .argsort()
                    .iloc[0]
                )
                idf = (
                    (df["Timestamp"] - faixa.iloc[i].to_numpy().max())
                    .abs()
                    .argsort()
                    .iloc[0]
                )
                dfn = pd.concat([dfn, df.iloc[range(idi, idf + 1)]], ignore_index=True)
            return dfn

        x = self.flow
        col_c = dd_c.flow_ext_df_name["marco"]  # = 'Marco'
        try:
            self.flow.flow_ext = grupa_m(x.flow_ext, col_c)
            self.flow.flow = ajuste(x.flow, self.flow.flow_ext)
            self.flow.flow_mean = grupa_m(x.flow_mean, col_c)
            self.flow.flow_std = grupa_m(x.flow_std, col_c)
            self.DTG = grupa_m(self.DTG, col_c)
            # TODO: Ta excluindo o que é inferior e superior
        except Exception:
            pass
        return self

    def select_DTG(self, X: List[str] = None, ID: List[int] = None) -> None:
        """Pega a distribuição do(s) local(is) desejado(s)

        Args:
            X (List[str]): Compare Linha principal.
            ID (List[int]): Identificação do local de extração.

            Opções para X: ['E_Well', 'E_ProdLine', 'E_ANM', 'E_FlowLine',
                    'E_Riser''E_Choke','E_Reservoir']
        """

        try:
            dfn = pd.DataFrame()
            xs = list([])
            dtgg = self.DTG.groupby("ID")
            self.compares = dict()
            for x in X:
                a, d = dd_b.X_LP[x]  # Pega o antes e depois
                self.compares[x] = [a, d]
                xs.append(a)
                xs.append(d)
                xs = list(set(xs))
            for x in xs:
                dfn = pd.concat([dfn, dtgg.get_group(x)])
        except Exception:
            print("Erro")
        self.DTG = dfn

    def get_DTG(self, teste: int = None, marco: int = None, ID: int = None):
        """Obtem a DTG específica de apenas uma análise

        Args:
            teste (int): Qual teste. Defaults to None.
            ID (int, optional): Identificacao do extrator. Defaults to None.

        Returns:
            _type_: _description_
        """
        col_b = dd_b.DTG_df_name["numero_escoamento"]  # = "N_escoam"
        dtgg = self.DTG.groupby(col_b).get_group(teste)
        dtgg_ID = dtgg.groupby("ID").get_group(ID)
        return dtgg_ID.groupby("Marco").get_group(marco)

    def calc_DP_GV(self) -> None:
        """Calcula a queda de pressao apenas na valvula

        Returns:
            _type_: _description_
        """
        # TODO: implementar isso
        if self.compares.keys() <= {"E_ANM"}:
            cols = dd_c.flow_circ_df_name
            # GV = dd_c.flow_circ_df_name['GV_ANM']
            # GV_position = self.flow.flow_mean[GV]
            DP = self.flow.flow_mean[cols["P1"]] - self.flow.flow_mean[cols["P2"]]
            self.flow.flow_mean["DP[bar]"] = DP

        elif self.compares.keys() <= {"E_Choke"}:
            cols = dd_c.flow_circ_df_name
            # GV = dd_c.flow_circ_df_name['GV_choke']
            # GV_position = self.flow.flow_mean[GV]
            # TODO: ajustar esse zero, pois não é zero
            DP = self.flow.flow_mean[cols["P4"]] - 0
            self.flow.flow_mean["DP[bar]"] = DP

        return self

    def get_prop(self, dir):
        """Obtem as propriedades da emulsao"""

        dir_mass = dir + "/massa_especifica.xlsx"
        dir_visc = dir + "/viscosidade.xlsx"
        # dir_tens = dir + "/tensaointerfacial.xlsx"

        emulsion = Fluid2(
            name="W/O",
            rho=pd.read_excel(dir_mass, sheet_name="emulsao"),
            mu=pd.read_excel(dir_visc, sheet_name="emulsao"),
            # sigma=pd.read_excel(dir_tens, sheet_name="emulsao"),
        )
        water = Fluid2(
            name="Saline Water",
            rho=pd.read_excel(dir_mass, sheet_name="agua"),
            mu=pd.read_excel(dir_visc, sheet_name="agua"),
        )
        oil = Fluid2(
            name="AW68",
            rho=pd.read_excel(dir_mass, sheet_name="oleo"),
            mu=pd.read_excel(dir_visc, sheet_name="oleo"),
        )
        self.tableprop = Prop(emulsion=emulsion, water=water, oil=oil)
        return self

    #
    #
    def preparaDados(self):
        """Seleciona em um só DataFrame, todas as informações necessárias
        para a solução do BP"""
        cols_ext = dd_c.flow_ext_df_name
        cols_int = dd_c.flow_circ_df_name
        Nemul = cols_ext["numero_emulsao"]
        Nesc = cols_ext["numero_escoamento"]
        ch2o = cols_ext["cH20"]
        marco = cols_ext["marco"]
        GV = cols_ext["GV_ANM"]
        T = cols_int["Tlinha"]
        dens = cols_int["densidade"]
        Vaz = cols_int["vazao"]

        self.dados = None
        self.dados = self.flow.flow_ext[[Nemul, Nesc, ch2o, marco, GV]]
        self.dados = pd.concat(
            [self.dados, self.flow.flow_mean[[T, dens, Vaz, "DP[bar]"]]], axis=1
        )

        # Propriedades do fluido
        # NOTE: interp está ignorando os dados fora da faixa
        def interpola_e_concatena(dados, df: Fluid2, C=None):
            """Interpola as propriedades

            Args:
                dados (_type_): dataframe com os dados
                df (Fluid2): fluido
                C (_type_, optional): Concentração. Defaults to None.

            Returns:
                _type_: Retorna os dados interpolados
            """
            if C is None:
                poly1 = polyfit(df.mu["T [°C]"], df.mu["mu [Pa.s]"], deg=1)
                mu = pd.DataFrame(
                    polyval(poly1, self.dados[T]),
                    columns=[df.name + " [Pa.s]"],
                )
                poly1 = polyfit(df.rho["T [°C]"], df.rho["rho [kg/m³]"], deg=3)
                rho = pd.DataFrame(
                    polyval(poly1, self.dados[T]),
                    columns=[df.name + " [kg/m³]"],
                )
            else:
                a, b = list(), list()
                for concentra, temp in zip(C, self.dados[T]):
                    a.append(
                        polyval(polyfit(df.mu["T [°C]"], df.mu[concentra], deg=1), temp)
                    )
                    b.append(
                        polyval(
                            polyfit(df.rho["T [°C]"], df.rho[concentra], deg=3), temp
                        )
                    )

                mu, rho = pd.DataFrame(a, columns=[df.name + " [Pa.s]"]), pd.DataFrame(
                    b, columns=[df.name + " [kg/m³]"]
                )
            dados = pd.concat([dados, mu, rho], axis=1)
            return dados

        self.dados = interpola_e_concatena(self.dados, self.tableprop.water)
        self.dados = interpola_e_concatena(self.dados, self.tableprop.oil)
        self.dados = interpola_e_concatena(
            self.dados, self.tableprop.emulsion, C=self.dados[ch2o]
        )

        self.dados["sigma [N/m]"] = 0.0045
        return self

    #
    def reduce_DTG(self, M: int, sl: slice):
        """Reduz as classes da DTG experimental pra corresponder à DTG numérica

        Args:
            M (_type_): Number of classes
            sl (_type_): slice for DataFrame rows
        """
        a = 0
        DTG = pd.DataFrame([])
        for t in self.teste:
            pass
        for t in self.marco:  # loop  testes
            for marco in self.marco[t]:  # loop   marcos de cada teste
                for comp in self.compares:  # loop extratores de cada marco
                    E_ant = self.compares[comp][0]  # Extrator antes
                    E_dep = self.compares[comp][1]  # Extrator depois

                    x = self.get_DTG(teste=t, marco=marco, ID=E_ant)[sl]
                    y = self.get_DTG(teste=t, marco=marco, ID=E_dep)[sl]

                    if len(x) and len(y) == M:
                        a = 1
                        DTG = pd.concat([DTG, x, y], axis=0, ignore_index=True)
                    else:
                        a = 0
                        print("Ja reduzido")
        if a:
            self.DTG = DTG


class DTG_experiment:
    def __init__(self, d32=0.32e-03, phi=0.117, U=2.72, theta=600.0):
        self.phi = phi
        self.theta = theta
        self.U = U
        self.d32 = d32

    def __repr__(self) -> str:
        phi = "{:.{}f}".format(100 * self.phi, 2)
        d32 = "{:.{}e}".format(1000 * self.d32, 2)
        return (
            f"concentration {phi}%  |  flow velocity {self.U} m/s  |  "
            + f"residence time {self.theta} s  |  SMD diameter {d32} mm"
        )


def get_location(folder, os_envi: str = "USERPROFILE"):
    home = Path.home()
    path = join(home, folder)
    path = join(environ["USERPROFILE"], folder)
    return path
