infolevel[Groebner]:=10;
Et_hat := [223967408523591254646-w_0, c*q*w_0*y_0-c*w_0*x_0*y_0+b*w_0+w_1, 254922694224434007520-z_0, -c*q*w_0*y_0+h*z_0+z_1, 162660288632176059331-k_0, k_1, -w_1-371557363636251837948136572846245215798350190780427097323004587532028871247800082, ((-x_1*y_0+y_1*(q-x_0))*w_0+w_1*y_0*(q-x_0))*c+b*w_1+w_2, beta*v_0*x_0+d*x_0-lm+x_1, -beta*v_0*x_0+a*y_0+y_1, 374995373184192105284752125379045631732493894738427417832194500901956238977278700-z_1, -q*(w_0*y_1+w_1*y_0)*c+h*z_1+z_2, -w_2+616406089539496933567575200339546400700960737765909159933208630896189410862732686023194388525279644303681211754594480965193275264476103643754, ((q*y_2-x_0*y_2-2*x_1*y_1-x_2*y_0)*w_0+(-2*w_1*x_1+w_2*(q-x_0))*y_0+2*y_1*w_1*(q-x_0))*c+b*w_2+w_3, beta*x_0*v_1+(beta*v_0+d)*x_1+x_2, -beta*v_0*x_1-beta*v_1*x_0+a*y_1+y_2, -k_0*y_0+u*v_0+v_1, -622109677272237935999102034631146555115641248728465499628402784148588177882603318710935485315380949729611545508294609855789896344516044634580-z_2, -q*(w_0*y_2+2*w_1*y_1+w_2*y_0)*c+h*z_2+z_3, -w_3-1022605132900407401382435538402252642484283327304933706983217483355779356610435226216516060103795745022828266611583216788512251843601424501689057037490049814124761449589255253251083083456159972809612368, ((q*y_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0)*w_0+(q*w_3-3*w_1*x_2-3*w_2*x_1-w_3*x_0)*y_0+(3*q*y_2-3*x_0*y_2-6*x_1*y_1)*w_1+3*y_1*w_2*(q-x_0))*c+b*w_3+w_4, (v_0*x_2+2*v_1*x_1+v_2*x_0)*beta+d*x_2+x_3, (-v_0*x_2-2*v_1*x_1-v_2*x_0)*beta+a*y_2+y_3, -k_0*y_1-k_1*y_0+u*v_1+v_2, 1032067268642457076210555033959571683131687081871074756978622926248938021549285098093380689736274603247474519026668545276347703251963706602138855315793890762225624749766060712762919610951053386985260960-z_3, -3*(y_2*w_1+1/3*y_3*w_0+w_2*y_1+1/3*w_3*y_0)*q*c+h*z_3+z_4, -w_4+1696481062696013622147401211483467003782265291247443491341259028826657743949748642133886416701545386145129743783008459262217847131342774597585458972953880884842769600542390097183739055944353138122324763665540287989251233704801240802160580855563125825090001815046, ((q*y_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0)*w_0+(q*w_4-4*w_1*x_3-6*w_2*x_2-4*w_3*x_1-w_4*x_0)*y_0+(4*q*y_3-4*x_0*y_3-12*x_1*y_2-12*x_2*y_1)*w_1+(-6*w_2*y_2-4*w_3*y_1)*x_0+(4*q*w_3-12*w_2*x_1)*y_1+6*y_2*q*w_2)*c+b*w_4+w_5, (v_0*x_3+3*v_1*x_2+3*v_2*x_1+v_3*x_0)*beta+d*x_3+x_4, (-v_0*x_3-3*v_1*x_2-3*v_2*x_1-v_3*x_0)*beta+a*y_3+y_4, -k_0*y_2-2*k_1*y_1-k_2*y_0+u*v_2+v_3, k_2, -1712178552941850007584126463596031768813933536369104531325945416395497192795049507975296446697760660357482666323203382252386852974823043123803123158351904512354607434985784249995208569928645844957031756332085352039532865320364809349202011886641898224394986797540-z_4, -4*q*(w_3*y_1+3/2*y_2*w_2+y_3*w_1+1/4*y_4*w_0+1/4*w_4*y_0)*c+h*z_4+z_5, -w_5-2814427488666333390517597550599129796697640955844878307189458035824586543913540705749117039953883490407298621950847650547506118583732185809731741244638228497297376418288569510087467844375292463153005893839220241441639388861849241261018048952076971961054887679593709853473508817392762673942342049705495947369021340084244912, ((q*y_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0)*w_0+(q*w_5-5*w_1*x_4-10*w_2*x_3-10*w_3*x_2-5*w_4*x_1-w_5*x_0)*y_0+(5*q*y_4-5*x_0*y_4-20*x_1*y_3-30*x_2*y_2-20*x_3*y_1)*w_1+(-10*w_2*y_3-10*w_3*y_2-5*w_4*y_1)*x_0+(5*q*w_4-30*w_2*x_2-20*w_3*x_1)*y_1+(10*w_2*y_3+10*w_3*y_2)*q-30*x_1*y_2*w_2)*c+b*w_5+w_6, (v_0*x_4+4*v_1*x_3+6*v_2*x_2+4*v_3*x_1+v_4*x_0)*beta+d*x_4+x_5, (-v_0*x_4-4*v_1*x_3-6*v_2*x_2-4*v_3*x_1-v_4*x_0)*beta+a*y_4+y_5, -k_0*y_3-3*k_1*y_2-3*k_2*y_1-k_3*y_0+u*v_3+v_4, k_3, 2840469304883688903673058097376076085467616154230446911648787529728574573384003999856918694225087522532578066525683969536578422612053177838653541081465202173609863266202724840704608379996031913577561766887816328675351014931830241639641727563999470617656690869853244504313524395405080397574741352817453934343385816879940980-z_5, -q*(w_0*y_5+5*w_1*y_4+10*w_2*y_3+10*w_3*y_2+5*w_4*y_1+w_5*y_0)*c+h*z_5+z_6, -w_6+4669077812382289007268257501248776221027459179049148874477498884223505966862603557216289653062101478258789724641790839365013967317927875820270475322877037407272127143057292587209308090684218867930234720529013526184371193850149463494140971267015692768745394849771522498790399925701050349583148246744316485608109937866810510944498092531918598000444211261046897777867126598738311870884, ((q*y_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0)*w_0+(q*w_6-6*w_1*x_5-15*w_2*x_4-20*w_3*x_3-15*w_4*x_2-6*w_5*x_1-w_6*x_0)*y_0+(6*q*y_5-6*x_0*y_5-30*x_1*y_4-60*x_2*y_3-60*x_3*y_2-30*x_4*y_1)*w_1+(-15*w_2*y_4-20*w_3*y_3-15*w_4*y_2-6*w_5*y_1)*x_0+(6*q*w_5-60*w_2*x_3-60*w_3*x_2-30*w_4*x_1)*y_1+(15*w_2*y_4+20*w_3*y_3+15*w_4*y_2)*q+(-60*x_1*y_3-90*x_2*y_2)*w_2-60*x_1*y_2*w_3)*c+b*w_6+w_7, (v_0*x_5+5*v_1*x_4+10*v_2*x_3+10*v_3*x_2+5*v_4*x_1+v_5*x_0)*beta+d*x_5+x_6, (-v_0*x_5-5*v_1*x_4-10*v_2*x_3-10*v_3*x_2-5*v_4*x_1-v_5*x_0)*beta+a*y_5+y_6, -k_0*y_4-4*k_1*y_3-6*k_2*y_2-4*k_3*y_1-k_4*y_0+u*v_4+v_5, k_4, -w_7-7745904879720683898539189396730195281728088268714593399757305518502085501653183450793781382129110000153553361947798529710894263906261470188793828229832061080610189453847732821699548951702904294321949984710886660851035678944883217601615730610571046755835066833248151567064454069839628511336841553736848550126441474210114818894917692101934055734046569902083021620214141510706554935638970826330884570504970347305258277959375582716086894721304138, -4712280654446701220551659336891300262878599901522794996998615162471427023161563336652025632816204811522460275485251936378895035017770494413370065008437054639243642783252361143442251745663698452981108719578708697544800725450203730442366681810274956922517704968334692098871730241945996457431525330987331845789630143537066020218730976150854317208809703700203953840265425502871324798020-z_6, -k_1, -k_2, -k_3, -k_4, 612808963620750075523655937147709763100254121029195352145001590373265572829792484850427015875728407-x_2, z_aux-1];
vars:=[w_7, z_6, y_6, x_6, w_6, z_5, y_5, x_5, w_5, v_5, z_4, y_4, x_4, w_4, v_4, k_4, z_3, y_3, x_3, w_3, v_3, k_3, z_2, y_2, x_2, w_2, v_2, k_2, z_1, y_1, x_1, w_1, v_1, k_1, z_0, y_0, x_0, w_0, v_0, k_0, z_aux, w_aux, a, b, beta, c, d, h, lm, q, u];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [k, x_2] 97
# 258, [1 = 16, v_4 = 1, a*y_0 = 1, b*w_1 = 1, -10*c*w_3*x_0*y_2 = 1, -w_0 = 1, -2*c*w_1*x_0*y_1 = 1, -4*c*q*w_1*y_3 = 1, -c*w_1*x_0*y_0 = 1, c*q*w_2*y_3 = 1, c*q*w_0*y_2 = 1, -5*c*w_1*x_0*y_4 = 1, -6*c*w_0*y_2 = 1, beta*v_0*x_5 = 1, -c*w_0*x_0*y_4 = 1, -beta*v_0*x_3 = 1, -c*w_0*x_5*y_0 = 1, c*q*w_1*y_2 = 1, u*v_4 = 1, -c*q*w_5*y_0 = 1, w_1 = 1, beta*v_3*x_1 = 1, -15*c*w_0*x_4*y_2 = 1, -y_0*w_0*c = 1, -c*w_0*x_1*y_0 = 1, -4*c*w_3*x_0*y_1 = 1, -2*c*w_1*x_1*y_0 = 1, -10*c*q*w_3*y_2 = 1, c*q*w_3*y_1 = 1, h*z_2 = 1, -6*c*w_0*x_5*y_1 = 1, -3*c*w_2*x_0*y_1 = 1, -10*v_3*beta = 1, b*w_6 = 1, -12*c*w_2*x_1*y_1 = 1, y_2 = 1, -30*c*w_4*x_1*y_1 = 1, c*q*w_5*y_1 = 1, -c*w_3*x_0*y_0 = 1, -15*c*w_4*y_0 = 1, c*q*w_2*y_4 = 1, -z_5 = 1, -15*c*w_4*x_0*y_2 = 1, -12*c*w_1*y_1 = 1, -3*c*w_2*x_1*y_0 = 1, -beta*v_5*x_0 = 1, b*w_2 = 1, -c*w_0*x_0*y_1 = 1, v_1 = 1, -z_3 = 1, -y_1 = 1, -w_6 = 1, -5*c*w_0*x_4*y_1 = 1, -c*q*w_0*y_4 = 1, -5*c*q*w_4*y_1 = 1, -4*c*w_0*x_1*y_3 = 1, beta*v_0*x_1 = 1, -20*c*w_3*x_1*y_1 = 1, -c*w_0*x_0*y_0 = 1, -c*w_0*x_6*y_0 = 1, b*w_3 = 1, z_2 = 1, y_3 = 1, -z_0 = 1, c*q*w_5*y_0 = 1, x_6 = 1, -10*c*q*w_2*y_3 = 1, z_5 = 1, -c*q*w_0*y_2 = 1, -c*w_0*x_3*y_0 = 1, v_3*beta = 1, -c*q*w_0*y_3 = 1, beta*v_1*x_1 = 1, -10*c*w_0*x_3*y_2 = 1, u*v_1 = 1, beta*v_1*x_4 = 1, -30*c*w_2*x_1*y_2 = 1, w_7 = 1, a*y_3 = 1, beta*x_0*v_1 = 1, -2*c*q*w_1*y_1 = 1, -w_2 = 1, z_4 = 1, -20*c*w_3*x_3*y_0 = 1, c*q*w_2*y_1 = 1, -c*q*w_4*y_0 = 1, beta*v_2*x_3 = 1, -60*c*w_3*y_1 = 1, -6*c*w_1*x_0*y_5 = 1, -w_5 = 1, -6*c*w_2*x_0*y_2 = 1, y_1 = 1, -6*c*q*w_2*y_2 = 1, c*q*w_2*y_2 = 1, -w_7 = 1, -3*c*w_1*x_0*y_2 = 1, -3*c*w_0*y_1 = 1, -10*c*w_2*x_3*y_0 = 1, -y_4 = 1, d*x_4 = 1, beta*v_2*x_1 = 1, c*q*w_6*y_0 = 1, v_2 = 1, -c*q*w_0*y_0 = 1, -60*c*w_1*y_3 = 1, c*q*w_1*y_5 = 1, -c*q*w_0*y_5 = 1, -c*q*w_1*y_0 = 1, a*y_1 = 1, -4*c*w_1*x_3*y_0 = 1, beta*v_2*x_0 = 1, -90*c*w_2*y_2 = 1, -c*q*w_0*y_1 = 1, c*q*w_0*y_1 = 1, -3*c*q*w_1*y_2 = 1, -30*c*w_1*x_1*y_4 = 1, v_2*beta = 1, x_3 = 1, -5*c*w_0*x_1*y_4 = 1, -w_3 = 1, -5*c*q*w_1*y_4 = 1, -c*w_2*x_0*y_0 = 1, -4*c*q*w_3*y_1 = 1, h*z_4 = 1, -3*c*w_1*y_0 = 1, -c*q*w_2*y_0 = 1, -6*c*w_1*x_1*y_1 = 1, -20*c*w_1*x_3*y_1 = 1, -3*c*q*w_2*y_1 = 1, -beta*v_0*x_4 = 1, c*q*w_0*y_4 = 1, -beta*v_0 = 1, -3*c*w_0*x_1*y_2 = 1, -20*c*w_3*x_0*y_3 = 1, -5*beta*v_1*x_4 = 1, -10*c*w_3*y_0 = 1, b*w_4 = 1, -z_6 = 1, c*q*w_4*y_2 = 1, -4*c*w_1*x_0*y_3 = 1, v_5 = 1, h*z_0 = 1, beta*v_0*x_3 = 1, -beta*v_0*x_1 = 1, c*q*w_0*y_6 = 1, -c*w_0*x_0*y_5 = 1, -c*w_0*x_0*y_6 = 1, c*q*w_4*y_0 = 1, d = 1, c*q*w_3*y_0 = 1, -12*c*w_1*x_1*y_2 = 1, beta*v_5*x_0 = 1, c*q*w_0*y_3 = 1, b*w_5 = 1, -c*w_0*x_0*y_3 = 1, beta*v_1*x_3 = 1, -5*c*w_1*x_4*y_0 = 1, -c*q*w_3*y_0 = 1, z_1 = 1, -y_3 = 1, -30*c*w_2*y_1 = 1, -5*c*w_4*x_1*y_0 = 1, -60*c*w_3*x_1*y_2 = 1, beta*v_0 = 1, beta*v_0*x_4 = 1, -w_4 = 1, -w_1 = 1, -6*v_2*beta = 1, c*q*w_3*y_2 = 1, -beta*x_0*v_1 = 1, beta*v_0*x_0 = 1, -30*c*w_1*x_4*y_1 = 1, x_1 = 1, beta*v_4*x_0 = 1, -15*c*w_2*x_4*y_0 = 1, c*q*w_1*y_1 = 1, c*q*w_0*y_5 = 1, -4*beta*v_1*x_3 = 1, -y_0 = 1, -y_2 = 1, -5*c*w_4*x_0*y_1 = 1, -z_1 = 1, -z_4 = 1, y_6 = 1, h*z_5 = 1, -10*c*w_2*x_0*y_3 = 1, c*q*w_0*y_0 = 1, x_4 = 1, beta*v_3*x_0 = 1, -beta*v_0*x_5 = 1, -15*c*w_0*y_4 = 1, -60*c*w_1*x_3*y_2 = 1, w_5 = 1, w_2 = 1, d*x_1 = 1, -6*c*w_5*x_0*y_1 = 1, -6*c*w_1*x_5*y_0 = 1, -beta*v_3*x_0 = 1, c*q*w_1*y_3 = 1, -6*c*w_2*y_0 = 1, h*z_1 = 1, -10*c*w_0*y_3 = 1, h*z_3 = 1, w_3 = 1, c*q*w_3*y_3 = 1, d*x_5 = 1, v_3 = 1, -beta*v_0*x_0 = 1, beta*v_1 = 1, -2*c*w_0*x_1*y_1 = 1, -20*c*w_0*x_3*y_3 = 1, x_5 = 1, b*w_0 = 1, -4*beta*v_3*x_1 = 1, -4*c*w_0*x_3*y_1 = 1, beta*v_4*x_1 = 1, -c*w_4*x_0*y_0 = 1, -5*beta*v_4*x_1 = 1, y_5 = 1, a*y_2 = 1, -3*beta*v_1 = 1, -3*beta*v_2*x_1 = 1, -20*c*w_1*x_1*y_3 = 1, y_4 = 1, w_6 = 1, -60*c*w_2*x_3*y_1 = 1, u*v_0 = 1, -10*beta*v_2*x_3 = 1, u*v_2 = 1, -c*w_0*x_0*y_2 = 1, a*y_5 = 1, w_4 = 1, a*y_4 = 1, c*q*w_4*y_1 = 1, -z_2 = 1, -60*c*w_2*x_1*y_3 = 1, -2*beta*v_1*x_1 = 1, -c*w_5*x_0*y_0 = 1, -4*c*w_3*x_1*y_0 = 1, -6*c*w_5*x_1*y_0 = 1, -beta*v_4*x_0 = 1, -c*w_0*x_4*y_0 = 1, d*x_0 = 1, z_3 = 1, d*x_3 = 1, u*v_3 = 1, -lm = 1, z_6 = 1, -beta*v_2*x_0 = 1, c*q*w_1*y_0 = 1, -c*w_6*x_0*y_0 = 1, -30*c*w_1*y_2 = 1, c*q*w_1*y_4 = 1, -6*c*w_0*x_1*y_5 = 1, c*q*w_2*y_0 = 1, z_aux = 1, -15*c*w_2*x_0*y_4 = 1, -1 = 1]
# 273, -5.446975752
# 1
# [x_2 = 25, k = 5]
# [x_2 = [4, 1, 4, 4, 3, 2, 3, 4, 4, 4, 3, 3, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 3, 3], k = [2, 2, 2, 2, 2]]
# [x_2 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], k = [1, 1, 1, 1, 1]]