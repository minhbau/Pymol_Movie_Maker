viewport 500, 500
mset 1x2200
mview store, 1, scene=F1
turn y, 120
mview store, 100, power = 1.0
turn y, 120
mview store, 200, power = 1.0
mview store, 300, scene=F1
mview store, 301, scene=F2
movie_fade cartoon_transparency, 301, 1.0, 302, 0.0
turn x, 120
mview store, 400, power = 1.0
turn x, 120
mview store, 500, power = 1.0
mview store, 600, scene=F2
mview store, 601, scene=F3
turn y, 120
mview store, 700, power = 1.0
turn y, 120
mview store, 800, power = 1.0
mview store, 900, scene=F3
movie_fade cartoon_transparency, 901, 0.0, 990, 1.0
mview store, 1000, scene=F5
turn x, 50
mview store, 1100, power = 1.0
turn x, 50
mview store, 1150, scene=F5
mview store, 1200, scene=F6
turn y, 50
mview store, 1300, power = 1.0
turn y, -100
mview store, 1450, power = 1.0
mview store, 1500, scene=F7
turn y, 60
mview store, 1700, power = 1.0
turn y, -120
mview store, 1800, power = 1.0
mview store, 1900, scene=F7
mview store, 1950, scene=F8
turn y, 60
mview store, 2000, power = 1.0
turn y, -120
mview store, 2100, scene=F8
mview store, 2200, scene=F8
mview reinterpolate
#set ray_trace_frames, 1