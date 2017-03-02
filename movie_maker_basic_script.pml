viewport 2000, 2000
mset 1x2000
mview store, 1, scene=F1
turn y, 120
mview store, 100, power = 1.0
turn y, 120
mview store, 200, power = 1.0
mview store, 300, scene=F1
turn x, 120
mview store, 400, power = 1.0
turn x, 120
mview store, 500, power = 1.0
mview store, 600, scene=F1
mview store, 601, scene=F2
turn x, 100
mview store, 700, power = 1.0
turn x, -100
mview store, 800, power = 1.0
mview store, 900, scene=F3
turn y, 120
mview store, 1000, power = 1.0
turn y, 120
mview store, 1100, power = 1.0
mview store, 1200, scene=F3
movie_fade cartoon_transparency, 1201, 0.0, 1290, 1.0
mview store, 1300, scene=F5
turn x, 50
mview store, 1400, power = 1.0
turn x, 50
mview store, 1500, scene=F5
mview store, 1600, scene=F6
mview store, 1700, scene=F7
turn y, 60
mview store, 1750, power = 1.0
turn y, -120
mview store, 1850, power = 1.0
mview store, 1900, scene=F7
mview store, 2000, scene=F7
mview reinterpolate
#set ray_trace_frames, 1