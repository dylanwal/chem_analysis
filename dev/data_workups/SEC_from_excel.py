import datetime

import numpy as np
import chem_analysis as ca

a = np.loadtxt(
    r"G:\Other computers\My Laptop\post_doc_2022\Data\polymerizations\DW2-10\DW2-10_RI.csv",
    delimiter=",",
               skiprows=2
)

x = a[:, 0]
data = a[:, 1::2]

time_ = [
'10:54:05',
'10:58:17',
'11:02:50',
'11:08:18',
'11:13:33',
'11:19:30',
'11:25:23',
'11:31:10',
'11:37:45',
'11:43:04',
'11:47:00',
'11:51:57',
'11:55:25',
'12:00:00',
'12:04:30',
'12:11:06',
'12:16:48',
'12:21:56',
'12:27:29',
'12:43:50',
'12:51:50',
'12:57:00',
'13:01:50',
'13:07:52',
'13:12:34',
'13:17:50',
'13:23:31',
'13:28:44',
'13:32:50',
'13:37:50',
'13:42:49',
'13:47:33',
'13:51:41',
'13:53:43',
'13:58:41',
'14:03:44',
'14:07:49',
'14:12:40',
'14:29:05',
'14:36:00',
'14:44:09',
'14:53:22',
'14:59:03',
'15:05:50',
'15:11:10',
'15:16:20',
'15:21:40',
'15:27:20',
'15:33:00',
'15:39:25',
'15:45:10',
'15:50:42',
'15:53:41',
'16:00:40',
'16:06:33',
'16:12:49',
'16:19:40',
'16:25:55',
'16:44:02',
'16:50:30',
'16:54:51',
'16:57:48',
'17:02:52',
'17:08:21',
'17:13:15',
'17:17:35',
'17:23:26',
'17:28:55',
'17:34:10',
'17:39:20',
'17:44:56',
'17:51:20',
'17:56:25',
'18:01:33',
'18:06:50',
'18:13:10',
'18:20:50',
'18:25:35',
'18:32:30',
'18:45:39',
'18:52:47',
'18:58:30',
'19:05:30',
    ]

time_ = np.array([datetime.datetime.fromisoformat(f"2023-12-07T{t}").timestamp() for t in time_])

d = ca.utils.general_math.pack_time_series(x, time_, data.T)

np.save(r"G:\Other computers\My Laptop\post_doc_2022\Data\polymerizations\DW2-10\DW2-10_RI", d)

print("hi")