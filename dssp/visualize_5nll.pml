load pdb5nll.ent.gz 
hide
show cartoon
set cartoon_fancy_helices,1
set cartoon_discrete_colors,1
set cartoon_highlight_color, grey60
set cartoon_dumbbell_length,1.0
set cartoon_rect_length,1.40000
set cartoon_loop_radius,0.3
set cartoon_smooth_loops=0
bg_color white
set_color helix = [0,1,0]
set_color strand = [1,0,0]
set_color coil = [1,1,0] 
select seg1, chain A and resi 1-6
select seg2, chain A and resi 7-9
select seg3, chain A and resi 10-27
select seg4, chain A and resi 27-28
select seg5, chain A and resi 29-34
select seg6, chain A and resi 35-39
select seg7, chain A and resi 40-46
select seg8, chain A and resi 47-48
select seg9, chain A and resi 48-55
select seg10, chain A and resi 56-63
select seg11, chain A and resi 64-77
select seg12, chain A and resi 78-79
select seg13, chain A and resi 80-89
select seg14, chain A and resi 89-92
select seg15, chain A and resi 93-107
select seg16, chain A and resi 107-180
select seg17, chain A and resi 109-118
select seg18, chain A and resi 119-123
select seg19, chain A and resi 124-138

color coil, seg2
color coil, seg4
color coil, seg6
color coil, seg8
color coil, seg10
color coil, seg12
color coil, seg14
color coil, seg16
color coil, seg18


color aquamarine, seg1
color hotpink, seg3
color deepolive, seg5
color dirtyviolet, seg7
color green, seg9
color firebrick, seg11
color pink, seg13
color deepsalmon, seg15
color purpleblue, seg17
color raspberry, seg19

#label the helices
#label /CA,"\316\261 "

#color strand, seg1
#color helix, seg3
#color strand, seg5
#color helix, seg7
#color strand, seg9
#color helix, seg11
#color strand, seg13
#color helix, seg15
#color strand, seg17
#color helix, seg19
set label_size, 30

label 17/CA,"\316\261A"
label 41/CA,"\316\261B"
label 64/CA,"\316\261C"
label 99/CA,"\316\261D"
label 125/CA,"\316\261E"

#label the strands
#label /CA,"\316\262 "
label 6/CA,"\316\2621"
label 34/CA,"\316\2622"
label 55/CA,"\316\2623"
label 88/CA,"\316\2624"
label 118/CA,"\316\2625"

alter A/1/, ss='S'
alter A/2/, ss='S'
alter A/3/, ss='S'
alter A/4/, ss='S'
alter A/5/, ss='S'
alter A/6/, ss='S'

alter A/10/, ss='H'
alter A/11/, ss='H'
alter A/12/, ss='H'
alter A/13/, ss='H'
alter A/14/, ss='H'
alter A/15/, ss='H'
alter A/16/, ss='H'
alter A/17/, ss='H'
alter A/18/, ss='H'
alter A/19/, ss='H'
alter A/20/, ss='H'
alter A/21/, ss='H'
alter A/22/, ss='H'
alter A/23/, ss='H'
alter A/24/, ss='H'
alter A/25/, ss='H'
alter A/26/, ss='H'
alter A/27/, ss='H'

alter A/29/, ss='S'
alter A/30/, ss='S'
alter A/31/, ss='S'
alter A/32/, ss='S'
alter A/33/, ss='S'

alter A/40/, ss='H'
alter A/41/, ss='H'
alter A/42/, ss='H'
alter A/43/, ss='H'
alter A/44/, ss='H'
alter A/45/, ss='H'
alter A/46/, ss='H'

alter A/49/, ss='S'
alter A/50/, ss='S'
alter A/51/, ss='S'
alter A/52/, ss='S'
alter A/53/, ss='S'
alter A/54/, ss='S'
alter A/55/, ss='S'

alter A/64/, ss='H'
alter A/65/, ss='H'
alter A/66/, ss='H'
alter A/67/, ss='H'
alter A/68/, ss='H'
alter A/69/, ss='H'
alter A/70/, ss='H'
alter A/71/, ss='H'
alter A/72/, ss='H'
alter A/73/, ss='H'
alter A/74/, ss='H'
alter A/75/, ss='H'
alter A/76/, ss='H'
alter A/77/, ss='H'

alter A/80/, ss='S'
alter A/81/, ss='S'
alter A/82/, ss='S'
alter A/83/, ss='S'
alter A/84/, ss='S'
alter A/85/, ss='S'
alter A/86/, ss='S'
alter A/87/, ss='S'
alter A/88/, ss='S'

alter A/93/, ss='H'
alter A/94/, ss='H'
alter A/95/, ss='H'
alter A/96/, ss='H'
alter A/97/, ss='H'
alter A/98/, ss='H'
alter A/99/, ss='H'
alter A/100/, ss='H'
alter A/101/, ss='H'
alter A/102/, ss='H'
alter A/103/, ss='H'
alter A/104/, ss='H'
alter A/105/, ss='H'
alter A/106/, ss='H'
alter A/107/, ss='H'

alter A/109/, ss='S'
alter A/110/, ss='S'
alter A/111/, ss='S'
alter A/112/, ss='S'
alter A/113/, ss='S'
alter A/114/, ss='S'
alter A/115/, ss='S'
alter A/116/, ss='S'
alter A/117/, ss='S'

alter A/124/, ss='H'
alter A/125/, ss='H'
alter A/126/, ss='H'
alter A/127/, ss='H'
alter A/128/, ss='H'
alter A/129/, ss='H'
alter A/130/, ss='H'
alter A/131/, ss='H'
alter A/132/, ss='H'
alter A/133/, ss='H'
alter A/134/, ss='H'
alter A/135/, ss='H'
alter A/136/, ss='H'
alter A/137/, ss='H'
alter A/138/, ss='H'

deselect
rebuild

