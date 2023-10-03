
class modelParameter:
    def __init__(self,
                 Cb: float = None,
                 Cc: float = None,
                 Ce: float = None,
                 varsigma: float = None) -> None:
        """Model parameters

        Args:
            Cb (float, optional): Breakage frequency model. Defaults to None.
            Cc (float, optional): coalescence frequency model. Defaults to None.
            Ce (float, optional): coalescence efficiency model. Defaults to None.
            varsigma (float, optional): mean number of daughter particles. Defaults to None.
        """
        self.Cb = Cb
        self.Cc = Cc
        self.Ce = Ce
        self.varsigma = varsigma
        self.modelparameters = {'Cb': Cb, 'Cc': Cc, 'Ce': Ce, 'varsigma': varsigma}
        #to get in a list: list(modelparameters.values())

    @property
    def parameters(self):
        #return [self._create_employee(**data) for data in self._employees]
        return list(self.modelparameters.values())

    def _add_parameters(self, *args):
        #self.modelparameters.update()
        pass

    # TODO: escrever função para adicionar parametros novos. 
    # TODO função para retornar modelparameters em lista