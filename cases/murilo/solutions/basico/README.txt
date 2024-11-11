murilo 29-04-2024:
	Todos os testes
	



C0_liao = [0.00481, 0.08, 2.8e-6, 1.83e9]

    (testes_all, "C&T Liao 5 ts", C0_liao, 5),  # 4
    (testes_all, "C&T Liao 10 ts", C0_liao, 10),
    (testes_all, "C&T Liao 20 ts", C0_liao, 20),
    (testes_all, "C&T Liao 40 ts", C0_liao, 40),
    (testes_all, "C&T Liao 50 ts", C0_liao, 50),
    (testes_all, "C&T Liao 60 ts", C0_liao, 60),
    (testes_all, "C&T Liao 100 ts", C0_liao, 100),
    (testes_all, "C&T Liao 150 ts", C0_liao, 150),  # 11
    (testes_all, "C&T Liao_0.1D", C0_liao, 100, 0.1),  # 12
    (testes_all, "C&T Liao_0.5D", C0_liao, 100, 0.5),
    (testes_all, "C&T Liao_1D", C0_liao, 100, 1),
    (testes_all, "C&T Liao_2D", C0_liao, 100, 2),
    (testes_all, "C&T Liao_2.4D", C0_liao, 100, 2.4),  # 16
    (testes_all, "C&T Liao_3D", C0_liao, 100, 3),
    (testes_all, "C&T Liao_4D", C0_liao, 100, 4),
    (testes_all, "C&T Liao_5D", C0_liao, 100, 5),
    (testes_all, "C&T Liao_6D", C0_liao, 100, 6),
    (testes_all, "C&T Liao_7D", C0_liao, 100, 7),  # 21
]
    "sol_C&T_liao_5ts": 4,
    "sol_C&T_liao_10ts": 5,
    "sol_C&T_liao_20ts": 6,
    "sol_C&T_liao_40ts": 7,
    "sol_C&T_liao_50ts": 8,
    "sol_C&T_liao_60ts": 9,
    "sol_C&T_liao_100ts": 10,
    "sol_C&T_liao_150ts": 11,
    "sol_C&T_liao_0.1D": 12,
    "sol_C&T_liao_0.5D": 13,
    "sol_C&T_liao_1D": 14,
    "sol_C&T_liao_2D": 15,
    "sol_C&T_liao_2.4D": 16,
    "sol_C&T_liao_3D": 17,
    "sol_C&T_liao_4D": 18,
    "sol_C&T_liao_5D": 19,
    "sol_C&T_liao_6D": 20,
    "sol_C&T_liao_7D": 21,

sol_Mitre_freqconst_S1_09-06-2024
    "P DDSD errado"
sol_Mitre_freqNconst_S1_09-06-2024
    "P DDSD errado"

sol_Mitre_freqconst_S1_10-06-2024
    testes_all,
    "Modelos do Mitre com frequencia de quebra constante",
    [1.0 / 1e2, 0.98 / 1e2],
    ["Cc", "Cb"],  # ordem das constantes
    {
        "breakup": "mitre_modified",
        "coalescence": "mitre_rigid_interface",
        "DDSD": "mitre",
        "varsigma": 26.0,  # S1: 26.0 ± 0.9. S3: 32.7 ± 16.8
    },
    100,

sol_Mitre_freqNconst_S1_10-06-2024
    testes_all,
    "Modelos do Mitre com frequencia não constante",
    [0.90 / 1e2, 1 / 1e2, 0.81 / 1e2],
    ["Cc", "Ce", "Cb"],  # ordem das constantes
    {
        "breakup": "mitre_modified",
        "coalescence": "mitre_partmobile_interface",
        "DDSD": "mitre",
        "varsigma": 31.2,
    },
    100,

sol_Mitre_freqconst_S3_10-06-2024
    testes_all,
    "Modelos do Mitre com frequencia de quebra constante",
    [1.88 / 1e2, 1.08 / 1e2],  # S3: [1.88 ± 0.06, 1.08 ± 0.07] . 10e-2
    ["Cc", "Cb"],  # ordem das constantes
    {
        "breakup": "mitre_modified",
        "coalescence": "mitre_rigid_interface",
        "DDSD": "mitre",
        "varsigma": 32.7,  # S1: 26.0 ± 0.9. S3: 32.7 ± 16.8
    },
    100