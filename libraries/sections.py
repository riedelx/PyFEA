import pickle
import numpy as np
import matplotlib.pyplot as plt
#import utils
#import math

class rss:
    def __init__(self, ID, mat, b, d):
        self.ID = ID
        self.mat = mat
        self.b = b
        self.d = d

    # def adaptic_print(self):
    #     # convert dictionary values of different dimensionality to a flatten list
    #     line = self.__dict__
    #     line = np.array(list(line.values()))
    #     return utils.str_joint(line)

class rccs:
    def __init__(self, ID, reinf_mat, unconf_mat, conf_mat, hc1, bc1, cover, links, reinf):
        self.ID = ID
        self.reinf_mat = reinf_mat
        self.unconf_mat = unconf_mat
        self.conf_mat = conf_mat
        self.hc1 = hc1
        self.hc2 = hc1 -2 * cover - links
        self.bc1 = bc1
        self.bc2 = bc1 -2 * cover - links
        # each layer of reinforcement given as: [no of bars, diameter, distance from the bottom fibre]
        # reinforcement given as a list [[layer_1], [layer_2],..., [layer_n]]
        # layers order starting from the furtherst layer from the bottom fibre
        self.reinf = [[int(math.pi*i[1]**2/4 * i[0]), i[2]] for i in reinf]
        self.reinf_basic = reinf # [no of bars, diameter, distance from the bottom fibre]
        self.area = hc1 * bc1
        self.cover = cover
        self.links = links # links diameter

    # def adaptic_print(self):
    #     # convert dictionary values of different dimensionality to a flatten list
    #     dict0 = self.__dict__
    #     dict1 = np.array(list(dict0.values())[0:8])
    #     dict2 = np.array(list(dict0.values())[8])
    #     dict2 = [val for sublist in dict2 for val in sublist]
    #     line = np.hstack((dict1,dict2))
    #     return utils.str_joint(line)

class rcts:
    def __init__(self, ID, reinf_mat, unconf_mat, conf_mat, Df, Dw, Bf, Bw, cover, links, reinf):
        self.ID = ID
        self.reinf_mat = reinf_mat
        self.unconf_mat = unconf_mat
        self.conf_mat = conf_mat
        self.Df = Df
        self.Dw = Dw
        self.df = Df - cover * 2 - links
        self.db = cover + Dw
        self.Bf = Bf
        self.Bw = Bw
        self.bf = Bf - 2 * cover - links
        self.bw = Bw - 2 * cover - links
        # each layer of reinforcement given as: [no of bars, diameter, distance from the bottom fibre]
        # reinforcement given as a list [[layer_1], [layer_2],..., [layer_n]]
        # layers order starting from the furtherst layer from the bottom fibre
        self.reinf = [[int(math.pi*i[1]**2/4 * i[0]), i[2]] for i in reinf]
        self.reinf_basic = reinf # [no of bars, diameter, distance from the bottom fibre]
        self.cover = cover
        self.links = links # links diameter
        area_w = Dw * Bw
        area_f = Df * Bf
        self.area = area_w + area_f
        self.centr_top = int((Df*area_f/2+area_w*(Df+Dw/2))/self.area)
        self.centr_bot = int(Df + Dw- self.centr_top)

    # def adaptic_print(self):
    #     # convert dictionary values of different dimensionality to a flatten list
    #     dict0 = self.__dict__
    #     dict1 = np.array(list(dict0.values())[0:12])
    #     dict2 = np.array(list(dict0.values())[12])
    #     dict2 = [val for sublist in dict2 for val in sublist]
    #     line = np.hstack((dict1,dict2))
    #     return utils.str_joint(line)
