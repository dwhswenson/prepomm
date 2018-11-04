"""
Tools to facilitate working with force fields in OpenMM
"""

class FFModels(object):
    def __init__(self, protein, water, **kwargs):
        self.protein = protein
        self.water = water
        self.kwargs = kwargs
        extras = [ff for ff in set(kwargs.values())
                  if ff not in [protein, water]]
        self.iter_list = [protein, water] + extras

    def __getitem__(self, item):
        try:
            return self.kwargs[item]
        except KeyError:
            return AttributeError()

    @property
    def water_model(self):
        """Water model name used by modeller"""
        return self.water.split('/')[-1].split('.')[0]

    def __iter__(self):
        return self.iter_list.__iter__()
