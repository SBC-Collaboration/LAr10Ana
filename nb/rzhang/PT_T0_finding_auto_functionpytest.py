import PT_T0_finding_auto_function as t0_pack
import PT_T0_finding_auto_function_v1 as t0_pack_v1
import numpy as np
# cp /exp/e961/data/SBC-25-daqdata/20251216_13.tar /exp/e961/app/users/runze/data/
# cd /exp/e961/app/users/runze/data/
# tar -xvf
test_path = f"/exp/e961/app/users/runze/data/20251211_6/"
save_path = test_path+f"fittings.txt"
# test_path = f"/exp/e961/app/users/runze/data/20251216_13/"
# t0_pack_v1.single_run_v1(test_path,2)
# t0_pack.single_run(test_path,2)

# v0 fitting function
# parameter_list = []
# for i in range(5,6,1):
#     try:
#         base_path = f"/exp/e961/app/users/runze/data/20251120_{i}/"
#         base_path = f"/exp/e961/app/users/runze/data/20251121_{i}/"
#         # base_path = f"/exp/e961/app/users/runze/data/20260112_{i}/"
#         # base_path = f"/exp/e961/app/users/runze/data/20251216_{i}/"
#         for j in range(0,99,1):
#             try:
#                 print(i,j, "new looop", base_path)
#                 paras = t0_pack.single_run(base_path,j)
#                 parameter_list.append(paras)
#             except Exception as e:
#                 print(e)
#                 continue
#     except:
#         continue
# print(parameter_list)

parameter_list = []
for i in range(5,6,1):
    try:
        base_path = f"/exp/e961/app/users/runze/data/20251120_{i}/"
        base_path = f"/exp/e961/app/users/runze/data/20251121_{i}/"
        # base_path = f"/exp/e961/app/users/runze/data/20260112_{i}/"
        # base_path = f"/exp/e961/app/users/runze/data/20251216_{i}/"
        for j in range(0,99,1):
            try:
                print(i,j, "new looop", base_path)
                paras = t0_pack_v1.single_run_v1(base_path,j)
                print("output", j, paras)
                parameter_list.append(paras)
            except Exception as e:
                print(e)
                continue
    except:
        continue
print(parameter_list)

parameter_array = np.array(parameter_list)
print(parameter_array[:,1])
print(parameter_array)

np.savetxt(save_path, parameter_array)

test_list=[(np.float64(1.4332005481360042e-05), np.float64(741.9772242140929)), (np.float64(1.482339818626237e-06), np.float64(610.5864684701319)), (np.float64(7.700009614105042e-06), np.float64(717.5004848709328)), (np.float64(2.386542092505951e-06), np.float64(665.3291558855778)), (np.float64(7.376689560893807e-06), np.float64(707.8019717291645)), (np.float64(3.223310034847078e-06), np.float64(688.2389873841931)), (np.float64(1.5247254629857806e-05), np.float64(725.2998066740577)), (np.float64(4.587644065881321e-06), np.float64(695.6386704191282)), (np.float64(3.630365356604633e-06), np.float64(663.0611814791816)), (np.float64(7.301147324369536e-06), np.float64(701.160871682847)), (np.float64(8.515245137912803e-07), np.float64(567.0153169713535)), (np.float64(5.382561696640237e-06), np.float64(696.6866140340255)), (np.float64(3.0788162209715117e-11), np.float64(837.8081812128087)), (np.float64(1.0501642687718843e-05), np.float64(697.596293430667)), (np.float64(5.764996671512247e-06), np.float64(714.5063132162112)), (np.float64(1.064633943018651e-06), np.float64(510.19374876162726)), (np.float64(0.0005730232671170386), np.float64(842.1838248995821)), (np.float64(2.917518958192357e-06), np.float64(660.7292918803344)), (np.float64(3.426953151626933e-06), np.float64(561.301967316283)), (np.float64(8.124380877376517e-06), np.float64(696.1124389339599))]

import numpy as np
import matplotlib.pyplot as plt
parameter_array = np.array(parameter_list)
fig, axes = plt.subplots(1,2)
axes[0].hist(parameter_array[:,0], bins=50, range=(0,5e-6))
axes[0].set_xlabel("a value")
axes[0].set_ylabel("counts")

axes[1].hist(parameter_array[:,1], bins=10)
axes[1].set_xlabel("t0 value/ms")
axes[1].set_ylabel("counts")
























