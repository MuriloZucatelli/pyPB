Dar uma geral nisso aqui pra por na dissertação de mestrado

opt_murilo_optglobal_shgo_erro_soma_29-05-2024:
	Testes selecionados
	função objetivo é a soma do erro de todos os testes
	iters=4,
	options={"ftol": 0.000001},
	bounds = [
        (1e-5 * 0.00481, 1e6 * 0.00481),
        (1e-6 * 0.08, 1e6 * 0.08),
        (1e-5 * 2.8e-6, 1e6 * 2.8e-6),
        (1e-6 * 1.83e9, 1e3 * 1.83e9),
    		]
	n = 2000
        method = "simplicial"


opt_murilo_optglobal_shgo_erro_soma_sobol_29-05-2024:
	Testes selecionados
	Função objetivo é a soma do erro de todos os testes
	iters=5,
	options={"ftol": 0.000001},
	bounds = [(1e-12, 1e12)] * 4
    		]
	n = 2000
        method = "sobol"
	Falhou


opt_murilo_optglobal_difevolution_erro_soma_29-05-2024
	Testes selecionados
	função objetivo é a soma do erro de todos os testes
	função evolução diferencial
	popsize = 400
	bounds = [
        (1e-5 * 0.00481, 1e6 * 0.00481),
        (1e-6 * 0.08, 1e6 * 0.08),
        (1e-5 * 2.8e-6, 1e6 * 2.8e-6),
        (1e-6 * 1.83e9, 1e3 * 1.83e9),
    		]
	 fun: 0.4866976561275103
	x: [ 2.691e+03  6.369e+04  2.428e-10  1.711e+12]
	nit: 44
	nfev: 72005
	time Pool:  125731.08595347404


opt_pymoo_CT_31-05-2024
	Testes selecionados
	função objetivo é a soma do erro de todos os testes
	função evolução diferencial
	popsize = 400
	# Constantes minimas
	Cl = array([1e-6 * 0.4, 1e-6 * 0.08, 1e-6 * 2.8e-6, 1e-6 * 1.83e9])
	# Constantes máximas
	Cu = array([1e6 * 0.4, 1e6 * 0.08, 1e6 * 2.8e-6, 1e4 * 1.83e9])
	Qual o resultado?????????????????
