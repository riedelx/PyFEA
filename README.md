# Introduction

PyFEA is a finite element analysis (FEA) software written in Python language. It is divided into separate analyses for different structural engineering calculations.

# Analyses

There are several available modules for different engineering calculations.

## 1. M-N interaction diagram 

Refer to file *mn-diagram.ipynb*.

This code generates moment interaction diagrams for reinforced concrete sections. 

### 1.1 Material models

This section defines how to incorporate material classes from _material.py_ module.

For each model, the parameters `plotting` and `title`  can be used to generate plots showing stress vs. strain diagrams.

#### 1.1.1 con1

```python
materials.con1(ID, fc1, length, epsilon_t2 = 0.001, fc2_factor = 0.1, ft_factor = 1, characteristic = True,plotting=True,title="stl1")
```

|   Parameter    |         Type         |                       Description                       |
| :------------: | :------------------: | :-----------------------------------------------------: |
|       ID       |         str          |                  name of the material                   |
|      fc1       |      float/int       |                peak compressive strength                |
|     length     |      float/int       |                     element length                      |
|   epsilon_t2   |        float         |                 ultimate tensile strain                 |
|   ft_factor    | float in range (0,1) |            tensile strength reduction factor            |
| characteristic |         bool         | characteristic strength (True) or mean strength (False) |

_Con1_ is a trilinear curve model in compression with an optional quadratic initial response . Tensile stage is given by a bilinear curve with softening.

![image](../assets/images/materials.con1_graph.png)

Initial compressive response is defined by the parameter $\alpha$, which is based on $E_{c0}$ and $E_{c1}$. Elastic initial modulus $E_{c0}$ is based on the parabolic curve. $E_{c1}$ is the secant modulus from the origin to the peak compressive stress. If $\alpha > 0$, a quadratic initial compressive response is implied.

After the peak compressive strength is achieved, softening stage takes place up to the failure. To avoid convergence issues, residual compressive strength $f_{c2}$ is maintained after failure. User can specify the _fc2_factor_ as the fraction of the peak strength. This is taken as 10% by default $f_{c2} = 0.1f_{c1}$.

The input parameters are concrete cylinder strength $f_{c1}$ and element length $h$. It is assumed that the input strength $f_{c1}$ is the characteristic compressive strength $f_{ck}$. If mean strength $f_{cm}$ is used, set the input parameter $characteristic$ to False, which affects the calculations of the fracture energy $G_f$. Element length $h$ is used to determine crack-band width.

Most of the other parameters are calculated according to <em>CEB-FIP Model Code 1990 (MC 1990)</em>, <em>CEB-FIP Model Code 2010 (MC 2010)</em> as well as <em>Rijkswaterstaat Technical Document: Guidelines for Nonlinear Finite Element Analysis of Concrete Structures (RTD 2010)</em>. These formulas are based on the uniaxial compressive cylinder strength.

To avoid overestimating the cracking moment, tensile strength $f_t$ can be reduced using tensile reduction factor *ft_factor*. Tension strain at failure $\varepsilon_{t2}$ needs to be defined by the user, taken as 0.001 by default.

|                  Parameter                   |                           Formula                            | Units |     Reference      |
| :------------------------------------------: | :----------------------------------------------------------: | :---: | :----------------: |
|        Compressive cylinder strength         |                  $$f_{c} = 0.85f_{c,cube}$$                  |  MPa  |         NA         |
| Characteristic compressive cylinder strength |                          $$f_{ck}$$                          |  MPa  |         NA         |
|          Mean compressive strength           |                   $$f_{cm} = f_{ck} + 8$$                    |  MPa  | MC 1990 Eq. 2.1-1  |
|          Peak compressive strength           |                          $$f_{c1}$$                          |  MPa  |         NA         |
|        Residual compressive strength         |                          $$f_{c2}$$                          |  MPa  |         NA         |
|               Tensile strength               | $$f_t= ft_{factor} \cdot 0.3f_{cm}^{2/3} \leq C50$$ $$ f_t= ft_{factor} \cdot 2.12ln(1+0.1f_{cm}) > C50$$ |  MPa  | MC 2010 Eq. 5.1-3a |
|               Fracture energy                |           $$G_f = 73\frac{ f_{cm}^{0.18}}{1000} $$           | N/mm  | MC 2010 Eq. 5.1-9  |
|         Initial compressive modulus          |           $$E_{c0} = 21500\cdot(f_{cm}/10)^{1/3}$$           |  MPa  | MC 2010 Eq. 5.1-21 |
|               Poisson's ratio                |                           $$0.2$$                            |   -   |  MC 2010 5.1.7.3   |
|         Compressive fracture energy          |                     $$G_{c} = 250G_{f}$$                     | N/mm  |   RTD 2010 p. 11   |
|     Compressive strain at peak strength      |      $$\varepsilon_{c1} = \frac{5}{3}\frac{f_c}{E_0}$$       |   -   |   RTD 2010 p. 21   |
|          Secant compressive modulus          |        $$E_{c1} =  \frac{f_{c1}}{\varepsilon_{c1}}$$         |  MPa  |         NA         |
|           Initial tensile modulus            |                     $$E_{t1} = E_{c0}$$                      |  MPa  |         NA         |
|          Compressive failure strain          | $$\varepsilon_{c2} = \varepsilon_{c1} + \frac{3}{2}\frac{G_c}{hf_c}$$ |   -   |   RTD 2010 p. 21   |
|       Tensile strain at peak strength        |          $$\varepsilon_{t1} = \frac{f_t}{E_{t1}}$$           |   -   |         NA         |
|            Tensile failure strain            |                     $$\varepsilon_{t2}$$                     |   -   |         NA         |
|     Initial compressive response factor      |          $$\alpha = \frac{E_{c0}-E_{c1}}{E_{c1}}$$           |   -   |   ADAPTIC manual   |

#### 1.1.2 stl1

```python
materials.stl1(ID, E1, fy, fu, epsilon_u,plotting=True,title="stl1",tension=True,compression=True)
```

| Parameter |   Type    |        Description        |
| :-------: | :-------: | :-----------------------: |
|    ID     |    str    |   name of the material    |
|    E1     | float/int | initial elastic stiffness |
|    fy     | float/int |      yield strength       |
|    fu     | float/int |     ultimate strength     |
| epsilon_u |   float   |      ultimate strain      |

_Stl1_ is a bilinear elasto-plastic model with kinematic strain hardening, used for a uniaxial modelling of mild steel. The curve can be defined for compression and/or tension by setting the parameters `compression`, `tension` to `True` or `False`.

<img src="../assets/images/materials.stl1.png" width="500" />

#### 1.1.3 EC2con

This is an elasto-perfectly-plastic material model, that works only in compression, defined in accordance to Eurocode 2 provisions.

`materials.EC2con(ID, fu, epsilon_u,plotting=True,title="EC2con",tension=True,compression=True)`

### 1.2 Sections

This chapter defines how to incorporate section classes from _sections.py_ module.

#### 1.2.1 rss

Rectangular solid section.

```python
sections.rss(ID, mat, b, d)
```

| Parameter | Type  |     Description     |
| :-------: | :---: | :-----------------: |
|    ID     |  str  | name of the section |
|    mat    |  str  |   material model    |
|     b     | float | section width [mm]  |
|     d     | float | section depth [mm]  |

<img src="../assets/images/sections.rss1.png" width="500" />

#### 1.2.2 rcts

Reinforced concrete T-section.

```python
sections.rcts(ID, reinf_mat, unconf_mat, conf_mat, Df, Dw, Bf, Bw, cover, links, reinf)
```

| Parameter  | Type  |           Description            |
| :--------: | :---: | :------------------------------: |
|     ID     |  str  |       name of the section        |
| reinf_mat  |  str  |   reinforcement material model   |
| unconf_mat |  str  | unconfined region material model |
|  conf_mat  |  str  |  confined region material model  |
|     Df     | float |        flange depth [mm]         |
|     Dw     | float |        web thickness [mm]        |
|     Bf     | float |        flange depth [mm]         |
|     Bw     | float |        web thickness [mm]        |
|   cover    | float |         cover depth [mm]         |
|   links    | float |       links diameter [mm]        |
|   reinf    | float |       reinforcement layers       |

Reinforcement layers _reinf_ have to be input as a nested list as follows:
```python
reinf = [layer_1,layer_2,...layer_n]
layer_k = [no_bars, dia, dist]
```
where the layer k has the following parameters:
* _no_bars_ - number of bars
* _dia_ - diamater of the bars [mm]
* _dist_ - distance from the bottom of the section [mm]

The layers have to be defined by _dist_ parameter in the descending order i.e. from the bars at the top of the section to the bottom.

<img src="../assets/images/sections.rcts1.png" width="400" />
<img src="../assets/images/sections.rcts2.png" width="250" />

