import sys, os
import subprocess
import time

script_dir = os.path.dirname(os.path.abspath(__file__))
print("Script directory:", script_dir)

args = sys.argv[1:]
if not args:
    print("Usage) python ./exp.py [category] [target]")
    print("\t 0[0~3]: SSD [18], Non-binary")
    print("\t 1[0~3]: SSD [18], Binary")
    print("\t 2[0~8]: PCG Beaver [22], Dual LPN, Binary")
    print("\t 3[0~4]: Committed VLOE [26,27], Primal, Binary")
    print("\t 4[0~1]: Ferret [31], Binary")
    print("\t 5[0~1]: Wolverine [30], Non-binary")
    print("\t 6[0]: TinyKey [15], Binary")
    print("\t 7[0]: MPCitH [9], Binary")
    print("\t 8[0~6]: PCGs [7,20], 80-bit security, Non-binary small weight")
    print("\t 9[0~6]: PCGs [7,20], 80-bit security, Non-binary Large weight")
    print("\t 10[0~6]: PCGs [7,20], 80-bit security, Binary small weight")
    print("\t 11[0~6]: PCGs [7,20], 80-bit security, Binary Large weight")
    print("\t 12[0~5]: PCGs [20], refined 128-bit security, Non-binary")
    print("\t 13[0~5]: PCGs [20], refined 128-bit security, Binary")
    sys.exit(1)

expCategory = [
    # 0) SSD [18], Non-binary
    [
    [2**10, 2**9, 107, 0],
    [2**12, 2**11, 90, 0],    
    [2**14, 2**13, 74, 0],
    [2**16, 2**15, 69, 0],
    ],
    # 1) SSD [18], Binary
    [
    [2**10, 2**9, 107, 1, 1.5],
    [2**12, 2**11, 90, 1],    
    [2**14, 2**13, 74, 1],
    [2**16, 2**15, 69, 1],
    ],
    # 2) PCG Beaver [22], Dual LPN, Binary
    [
    [2**20, 2**18, 47, 1],
    [2**21, 2**18, 31, 1],
    [2**21, 2**19, 45, 1],
    [2**22, 2**19, 30, 1],
    [2**22, 2**20, 44, 1],
    [2**23, 2**21, 42, 1],
    [2**24, 2**22, 41, 1],
    [2**25, 2**23, 39, 1],
    [2**26, 2**24, 37, 1],
    ],
    # 3) Committed VLOE [26,27], Primal, Binary
    [
    [ 1456, 760, 91, 1, 1.5],
    [ 1472, 762, 92, 1, 1.5],
    [ 1488, 765, 93, 1, 1.5],
    [ 1272, 433, 159, 1, 1.5],
    [ 1296, 434, 162, 1, 1.5],   
    ],
    # 4) Ferret [31], Binary
    [
    [609728, 36288, 1269, 1],
    [10805248, 589760, 1319, 1]
    ],
    # 5) Wolverine [30], Non-binary
    [
    [642048, 19870, 2508, 0],
    [10805248, 589760, 1319, 0]
    ],
    # 6) TinyKey [15], Binary
    [
    [1280,860,80,1,2.5],
    ],
    # 7) MPCitH [9], Binary
    [
    [1736, 776, 217, 1, 1.5]
    ],
    # 8) PCGs [7,20], 80-bit security, Non-binary small weight
    [
    [2**10, 652, 57, 0],
    [2**12, 1589, 98, 0],
    [2**14, 3482, 198, 0],
    [2**16, 7391, 389, 0],
    [2**18, 15336, 760, 0],
    [2**20, 32771, 1419, 0],
    [2**22, 64770, 2735, 0]
    ],
    # 9) PCGs [7,20], 80-bit security, Non-binary Large weight
    [
    [2**10, 652, 106, 0, 2],
    [2**12, 1589, 172, 0],
    [2**14, 3482, 338, 0],
    [2**16, 7391, 667, 0],
    [2**18, 15336, 1312, 0],
    [2**20, 32771, 2467, 0],
    [2**22, 64770, 4788, 0]
    ],
    # 10) PCGs [7,20], 80-bit security, Binary small weight
    [
    [2**10, 652, 57, 1, 1.5],
    [2**12, 1589, 98, 1],
    [2**14, 3482, 198, 1],
    [2**16, 7391, 389, 1],
    [2**18, 15336, 760, 1],
    [2**20, 32771, 1419, 1],
    [2**22, 64770, 2735, 1]
    ],
    # 11) PCGs [7,20], 80-bit security, Binary Large weight
    [
    [2**10, 652, 106, 1, 2.5],
    [2**12, 1589, 172, 1],
    [2**14, 3482, 338, 1],
    [2**16, 7391, 667, 1],
    [2**18, 15336, 1312, 1],
    [2**20, 32771, 2467, 1],
    [2**22, 64770, 4788, 1]
    ],
    # 12) PCGs [20], refined 128-bit security, Non-binary
    [
    [2**12,1377,172,0],
    [2**14,2909,338,0],
    [2**16,6091,667,0],
    [2**18,14796,1312,0],    
    [2**20,30978,2467,0],
    [2**22,75396,4788,0],
    ],
    # 13) PCGs [20], refined 128-bit security, Binary
    [
    [2**12,1657,172,1],
    [2**14,3655,338,1],
    [2**16,7560,667,1],
    [2**18,15996,1312,1],    
    [2**20,33354,2467,1],
    [2**22,80074,4788,1],
    ],
]

nameTag = [
    "0_SSD_Non-binary",
    "1_SSD_Binary",
    "2_PCG_Beaver",
    "3_Committed_VLOE",
    "4_Ferret",
    "5_Wolverine",
    "6_TinyKey",
    "7_MPCitH"
    "8_PCGs_80-bit_security_Non-binary_small_weight",
    "9_PCGs_80-bit_security_Non-binary_Large_weight",
    "10_PCGs_80-bit_security_Binary_small_weight",
    "11_PCGs_80-bit_security_Binary_Large_weight",
    "12_PCGs_128-bit_security_Non-binary",
    "13_PCGs_128-bit_security_Binary",
]
if len(args) == 1:
    exps = expCategory[int(args[0])]
    output_file = f"output_" + nameTag[int(args[0])] + ".txt"  # File
else:
    exps = expCategory[int(args[0])][int(args[1]):int(args[1])+1]
    output_file = f"output_" + nameTag[int(args[0])] + ".txt"  # File

env = os.environ.copy()
env["OMP_NUM_THREADS"] = "16"

with open(output_file, "w", encoding="utf-8") as f:
    for exp in exps:
        if len(exp) == 5:
            command = [script_dir + "/rsd", f"{exp[0]}", f"{exp[1]}", f"{exp[2]}", f"{exp[3]}", f"{exp[4]}"]
        else:
            command = [script_dir + "/rsd", f"{exp[0]}", f"{exp[1]}", f"{exp[2]}", f"{exp[3]}"]
                       
        f.write(f"================== ruinng RSD({exp[0]},{exp[1]},{exp[2]}) ==================\n")
        f.flush()
        print(f"=== ruinng RSD({exp[0]:7d},{exp[1]:7d},{exp[2]:7d}) ===")

        sp = time.perf_counter()
        result = subprocess.run(
            command,
            stdout=f,
            stderr=subprocess.STDOUT,
            text=True,
            env=env
        )
        ep = time.perf_counter()
        f.write(f"time : {ep-sp:.4f}s\n")
        f.write(f"return code: {result.returncode}\n\n")
        f.flush()
print("End")