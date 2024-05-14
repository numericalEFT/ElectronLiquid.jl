
import re
import matplotlib.pyplot as plt
import numpy as np

filename = "ver4.dat"

direct = []
exchange = []

direct_e = []
exchange_e = []

current_section = None

with open(filename, 'r') as file:
    for line in file:
        line = line.strip()
        if not line:  # Ignore empty lines
            continue
        if line.startswith("# Direct ver4"):
            current_section = "direct"
        elif line.startswith("# Exchange ver4"):
            current_section = "exchange"
        elif current_section:
            # line = re.sub(r'±\s+\S+', '', line)
            # print(re.split(r'\s+', line))
            angle, t, e1, te, u, e2, ue, s, e3, se, total, e4, tote = re.split(r'\s+', line)
            angle = float(angle)
            t = float(t.strip())
            u = float(u.strip())
            s = float(s.strip())
            total = float(total.strip())
            te = float(te.strip())
            ue = float(ue.strip())
            se = float(se.strip())
            tote = float(tote.strip())

            if current_section == "direct":
                direct.append((angle, t, u, s, total))
                direct_e.append((angle, te, ue, se, tote))
            elif current_section == "exchange":
                exchange.append((angle, t, u, s, total))
                exchange_e.append((angle, te, ue, se, tote))


direct = np.array(direct)
exchange = np.array(exchange)
direct_e = np.array(direct_e)
exchange_e = np.array(exchange_e)

print("Direct ver4:", direct)
print("Exchange ver4:", exchange)

# plt.errorbar(direct[:, 0], direct[:, 4], yerr=direct_e[:, 4], fmt='o', markersize=4, label='Direct Total', color='blue')
# plt.errorbar(exchange[:, 0], exchange[:, 4], yerr=exchange_e[:, 4], fmt='o', markersize=4, label='Exchange Total', color='red')
# plt.show()


s = direct[:, 4] + exchange[:, 4]/2
a =  exchange[:, 4]/2
# uu = s+a
# ud = s-a
angle = direct[:, 0]
# angle = np.sin(angle/2)**2

plt.errorbar(angle, s, yerr=direct_e[:, 4], fmt='o', markersize=4, label='UU Total', color='blue')
plt.errorbar(angle, a, yerr=exchange_e[:, 4], fmt='o', markersize=4, label='UD Total', color='red')
plt.legend()
plt.show()




