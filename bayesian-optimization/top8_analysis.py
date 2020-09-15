import numpy as np
from scipy import stats

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), stats.sem(a)
    h = se * stats.t.ppf((1 + confidence) / 2., n-1)
    return m, h

data_log = np.loadtxt('output_raw.log')

sorted_data_log = data_log[data_log[:,8].argsort()]
reversed_sorted_data_log = sorted_data_log[::-1]
top_8_data_log = reversed_sorted_data_log[:8]

ex1, p_ex1 = top_8_data_log[:, 0], top_8_data_log[:, 1]
ex2, p_ex2 = top_8_data_log[:, 2], top_8_data_log[:, 3]
ex3, p_ex3 = top_8_data_log[:, 4], top_8_data_log[:, 5]
background = top_8_data_log[:, 6]
mse = top_8_data_log[:, 7]
RunNo = top_8_data_log[:, 9]

# Sort so ex1 > ex2 > ex3 (In terms of q-value)
ex1_ = []
p_ex1_ = []

ex2_ = []
p_ex2_ = []

ex3_ = []
p_ex3_ = []

for i in range(len(ex1)):
    states = [(ex1[i], p_ex1[i]), (ex2[i], p_ex2[i]), (ex3[i], p_ex3[i])]
    sorted_states = sorted(states)
    ex1_.append(sorted_states[0][0])
    p_ex1_.append(sorted_states[0][1])
    ex2_.append(sorted_states[1][0])
    p_ex2_.append(sorted_states[1][1])
    ex3_.append(sorted_states[2][0])
    p_ex3_.append(sorted_states[2][1])
    print(RunNo[i],mse[i],ex1_[i],p_ex1_[i],ex2_[i],p_ex2_[i],ex3_[i],p_ex3_[i],background[i])


background_m, background_h = mean_confidence_interval(background)
print("Background: {:5.3f} +/- {:5.3f}".format(background_m, background_h))
#±

ex1_m, ex1_h = mean_confidence_interval(ex1_)
print("ex1: {:5.3f} +/- {:5.3f}".format(ex1_m, ex1_h))

ex2_m, ex2_h = mean_confidence_interval(ex2_)
print("ex2: {:5.3f} +/- {:5.3f}".format(ex2_m, ex2_h))

ex3_m, ex3_h = mean_confidence_interval(ex3_)
print("ex3: {:5.3f} +/- {:5.3f}".format(ex3_m, ex3_h))

p_ex1_m, p_ex1_h = mean_confidence_interval(p_ex1_)
print("p_ex1: {:5.3f} +/- {:5.3f}".format(p_ex1_m, p_ex1_h))

p_ex2_m, p_ex2_h = mean_confidence_interval(p_ex2_)
print("p_ex2: {:5.3f} +/- {:5.3f}".format(p_ex2_m, p_ex2_h))

p_ex3_m, p_ex3_h = mean_confidence_interval(p_ex3_)
print("p_ex3: {:5.3f} +/- {:5.3f}".format(p_ex3_m, p_ex3_h))
