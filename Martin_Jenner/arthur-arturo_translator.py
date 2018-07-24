import cPickle
import numpy.core.multiarray
import matplotlib.pyplot as plt
import numpy as np
import collections
plt.close(1)
plt.close(2)
plt.close(3)
docpath = r"C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\NeoNarval\Martin_Jenner\test_ThAr\Bmatrix_data_sheet.txt"
dic = collections.OrderedDict()
with open(docpath, 'r') as file:
    content = file.readlines()
content = [line.strip() for line in content]
for line in content:
    param = line.split(" : ")
    if len(param) == 2:
        nom_param = param[0]
        value_param = param[1]
        dic[nom_param] = value_param
pix = float(dic["pix/arturo"])
lambd_ini = int(dic["initial index"])*pix
lambd_end = int(dic["final index"])*pix
order = dic["order"]
lane = dic["lane"]
path = r"C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\TEMP_\ThAr_based_spec_OR"+order+"_LA"+lane+".p"
print(path)

order_ref = dic["order ref"]
order_curv = dic["order curv"]
order_ref_data = cPickle.load(open(order_ref, 'r'))
order_curv_data = cPickle.load(open(order_curv, 'r'))

#----fichier type lanes_thickness--------
thic_data = []
B = order_ref_data['Beams']
I = order_ref_data['Orders']
for j in range(40):
    for i in range(2):
        if (2*j+i)%2 == 0:
            thic_data.append(int(B[j]))
        else:
            thic_data.append(int(I[j+1]-I[j]-B[j]))
            
cPickle.dump(thic_data, open(r'C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\NeoNarval\Martin_Jenner\arturo_lanes_thic.p', 'w'))
#----------------------------------------

#----fichier type lanes_position---------

map = np.zeros([4612,80])

shift = order_curv_data['map'].T
newshift = np.zeros([4612,80])
P = []
for j in range(40):
    for i in range(2):
        if (2*j+i)%2 == 0:
            P.append(I[j])
            newshift[:,2*j+i]=shift[:,j]
        else:
            P.append(I[j]+B[j])
            newshift[:,2*j+i]=shift[:,j]
map[2306,:] = P
for k in range(1,2306):
    map[2306+k,:] = P + newshift[2306+k,:]
    map[2306-k,:] = P + newshift[2306-k,:]

cPickle.dump(map, open(r'C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\NeoNarval\Martin_Jenner\arturo_lanes_pos.p', 'w'))
#----------------------------------------

#----fichier type poly_envelope----------
C = order_curv_data['Coefs']
poly_data = []
for j in range(39):
    for i in range(2):
        poly_data.append(np.poly1d(np.array(C[j]) + np.array([0,0,0,I[j]])))

cPickle.dump(poly_data, open(r'C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\NeoNarval\Martin_Jenner\arturo_poly_env.p', 'w'))
#----------------------------------------
























