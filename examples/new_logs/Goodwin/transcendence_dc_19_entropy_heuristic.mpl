infolevel[Groebner]:=10;
Et_hat := [3912304112999020354-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, 746715582683367512-gama_0, gama_1, -x1_1-30279812051594593527571393129418542111497234109986194399/2863022383777748760, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama_0*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+30237858788335825459088492773507940259219757263269502699656852899889712368366919116802920017141451823019246589/1057612992444663784681488214706092918777015965932935200, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_1-gama_0*x2_1-gama_1*x2_0)*x4_0-x4_1*(-delta*x3_0+gama_0*x2_0))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama_0*x2_0+x3_1, -x1_3-15097981826065988168566958485245833363468589263593933395283300987795222838812995953093908539696245974694540604889023715695930414795544818378336544066797162079552319/195343432892033423946854714170238895106306238577542896913147739154239090146076448688952000, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0)*x4_0+(2*delta*x3_1-2*gama_0*x2_1-2*gama_1*x2_0)*x4_1+x4_2*(delta*x3_0-gama_0*x2_0))*sgm+2*x3_1*x4_2+x3_0*x4_3+x3_2*x4_1, gama_2, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama_0*x2_1-gama_1*x2_0+x3_2, -x1_4-70679414757638856857315547122997489309595970807533110945020524821609462573365417189408754058530123627394100525680388813839851041742770047718477354412915992559946827029712185068511137456936665735301537275245768759369/2672619207848142107756947426085631186747526245932866392132165719734558248530911875008485332814908415168340524720246391520000, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0)*x4_0+(3*delta*x3_2-3*gama_0*x2_2-6*gama_1*x2_1-3*gama_2*x2_0)*x4_1+(-x2_0*x4_3-3*x2_1*x4_2)*gama_0+(x3_0*x4_3+3*x3_1*x4_2)*delta-3*x4_2*gama_1*x2_0)*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, gama_3, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0+x3_3, -x1_5-281243124771611200267571225815260511239944289340385272844174274381792026073548913922231949580637076810250454991180044377060040979258374895527064359745354700222905527001633674997371737342264103807664668628371491722626943486418857719093680558997051238205147138311332214236934004017806372640593/4442747518643056502182860993439889241884708823516259209347745249674746567100358504449890002333351018120880784479822822450773058269705519643452210298009876800000, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0)*x4_0+(4*delta*x3_3-4*gama_0*x2_3-12*gama_1*x2_2-12*gama_2*x2_1-4*gama_3*x2_0)*x4_1+(-x2_0*x4_4-4*x2_1*x4_3-6*x2_2*x4_2)*gama_0+(x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2)*delta+(-4*gama_1*x4_3-6*gama_2*x4_2)*x2_0-12*x4_2*gama_1*x2_1)*sgm+4*x3_3*x4_2+x3_4*x4_1+6*x3_2*x4_3+4*x3_1*x4_4+x3_0*x4_5, gama_4, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0+x3_4, -x1_6-241091338228430857671060885086519055166184409350403301364383755028388078266336222102806698917192215814169664295754199047522794986493102291261070264463991999748613509825573047797709401210078082181501480634112170472732963365458589051405186538441677178619946751231939488757399984042773293133733273775891441138960245674472606199382755260006542320251735409476713807043/164117036754293844269170621090763024493311375123794216770137339573687773610723715330707121534339148630535983625189051731282263747750977240338287368819982973346256215107548473583988596033033600000, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0)*x4_0+(5*delta*x3_4-5*gama_0*x2_4-20*gama_1*x2_3-30*gama_2*x2_2-20*gama_3*x2_1-5*gama_4*x2_0)*x4_1+(-x2_0*x4_5-5*x2_1*x4_4-10*x2_2*x4_3-10*x2_3*x4_2)*gama_0+(x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2)*delta+(-5*gama_1*x4_4-10*gama_2*x4_3-10*gama_3*x4_2)*x2_0+(-30*gama_1*x2_2-30*gama_2*x2_1)*x4_2-20*x4_3*gama_1*x2_1)*sgm+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+5*x3_1*x4_5+x3_0*x4_6, gama_5, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0+x3_5, -x1_7+142660373865189530707606541454663957314316125829783727459881349421098073419824576574155119829485417234031616507328456183253307161114934694127045123849884778707847725422442652220961147722588375500504357656667017945742764266179952707573553292601505293713148332948364468937811857594338860216606049262088171816052083584950315240655582969722479350698702616001518425874177860321849453360552667543696642921607941124131021932715481535804064589/2526064636802796650887165341412252201753831744595654257120123383138634786804040705709287233985258582830248889956316890409517948111597670301822958220246353088390252365011637468514616358841804640228912179002967483313987991928000000, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0)*x4_0+(6*delta*x3_5-6*gama_0*x2_5-30*gama_1*x2_4-60*gama_2*x2_3-60*gama_3*x2_2-30*gama_4*x2_1-6*gama_5*x2_0)*x4_1+(-x2_0*x4_6-6*x2_1*x4_5-15*x2_2*x4_4-20*x2_3*x4_3-15*x2_4*x4_2)*gama_0+(x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2)*delta+(-6*gama_1*x4_5-15*gama_2*x4_4-20*gama_3*x4_3-15*gama_4*x4_2)*x2_0+(-60*gama_1*x2_3-90*gama_2*x2_2-60*gama_3*x2_1)*x4_2+(-30*x2_1*x4_4-60*x2_2*x4_3)*gama_1-60*gama_2*x2_1*x4_3)*sgm+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1+x3_0*x4_7, gama_6, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0+x3_6, -x1_8+10280693357440931046460774432978350223164445036023878590999820015121728327218104397333488311248156814153300370063369935395051501354515293302261126086353619386925442096500795187544560152315170244512687480006810958893097652219875036054891881550930228659035943749611374519123176318416382635081165039311909404428083940875287869020594068400101755225240257168495503930507915216431202721356813940436003013100884822570058285841515073597014221430802017748807658988558313016085273622920177119727843597828656017283838599/1866278653478414379119788848817650010125240073414940926862198219987176358856497213057564996273925746311898239895089165469198548280650813129473274881400242476330993937960400194733690551813810425609139726884318273457493109067409921066180809625490778518863053120000000, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((delta*x3_7-gama_0*x2_7-7*gama_1*x2_6-21*gama_2*x2_5-35*gama_3*x2_4-35*gama_4*x2_3-21*gama_5*x2_2-7*gama_6*x2_1-gama_7*x2_0)*x4_0+(7*delta*x3_6-7*gama_0*x2_6-42*gama_1*x2_5-105*gama_2*x2_4-140*gama_3*x2_3-105*gama_4*x2_2-42*gama_5*x2_1-7*gama_6*x2_0)*x4_1+(-x2_0*x4_7-7*x2_1*x4_6-21*x2_2*x4_5-35*x2_3*x4_4-35*x2_4*x4_3-21*x2_5*x4_2)*gama_0+(x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2)*delta+(-7*gama_1*x4_6-21*gama_2*x4_5-35*gama_3*x4_4-35*gama_4*x4_3-21*gama_5*x4_2)*x2_0+(-105*gama_1*x2_4-210*gama_2*x2_3-210*gama_3*x2_2-105*gama_4*x2_1)*x4_2+(-42*x2_1*x4_5-105*x2_2*x4_4-140*x2_3*x4_3)*gama_1+(-105*gama_2*x4_4-140*gama_3*x4_3)*x2_1-210*gama_2*x2_2*x4_3)*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, gama_7, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0+x3_7, -x1_9+10848987737445746196663886156398937432655434952864502826391320394207675183818285403039733686242448729606540629993348042489516009252154466916932555595875067124197762793445298022744696017783546231301081234622956555674529076866347814202954902015447947821736913312249942723920519686199895666099558366135297186310150316282972741576902486213694947056624915212561287211383428596353584356643592824952371764157733048243360923350905284419028868746411463239126435470215298383759565684783361006433777886661013392402643047713966553245431161299013817554701664756912981351470684703007363178681717/172352874590215405090563647619276032597633261909320322956002238164961086138097489133695570495342711303739480007650032083472040116636838725946039660469148015525554571330173003041750586802919061023142145182374411693927926583059462542071254085665321778324435202689628411504143892652969266356145600000000, -gama_1, -gama_2, -gama_3, -gama_4, -gama_5, -gama_6, -gama_7, -8023727011106490690987333399344306890-x2_1, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, gama_7, x4_7, x3_7, x2_7, x1_7, gama_6, x4_6, x3_6, x2_6, x1_6, gama_5, x4_5, x3_5, x2_5, x1_5, gama_4, x4_4, x3_4, x2_4, x1_4, gama_3, x4_3, x3_3, x2_3, x1_3, gama_2, x4_2, x3_2, x2_2, x1_2, gama_1, x4_1, x3_1, x2_1, x1_1, gama_0, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [gama, x2_1] 191
# 274, [1 = 7, -1 = 3, x3_1*x4_5 = 1, delta*x3_2 = 1, x2_2 = 1, x3_1 = 1, x2_5 = 1, x1_3*x4_1 = 1, -21*sgm*x2_5*x4_2 = 1, b*x1_0*x4_5 = 1, x1_2*x4_7 = 1, x3_4*x4_3 = 1, -3*sgm*x2_2*x4_1 = 1, x1_1*x4_0 = 1, x1_1*x4_5 = 1, delta*x3_0 = 1, delta*sgm*x3_0*x4_1 = 1, b*x1_5*x4_0 = 1, x1_7*c = 1, delta*sgm*x3_0*x4_0 = 1, -x1_5 = 1, b*x1_2*x4_3 = 1, x2_7 = 1, x1_2*x4_5 = 1, -sgm*x2_5*x4_0 = 1, x2_6 = 1, -x1_9 = 1, beta*x2_6 = 1, delta*sgm*x3_3*x4_4 = 1, delta*sgm*x3_2*x4_4 = 1, x3_0*x4_3 = 1, z_aux*x3_0*x4_0 = 1, -5*sgm*x4_4 = 1, x3_4 = 1, delta*sgm*x3_4*x4_1 = 1, x3_1*x4_4 = 1, b*c*x1_7 = 1, b*x1_1*x4_3 = 1, -3*sgm*x4_2 = 1, x3_5*x4_2 = 1, -6*sgm*x2_2*x4_2 = 1, delta*x3_6 = 1, x1_3*x4_6 = 1, b*x1_6*x4_0 = 1, delta*sgm*x3_3*x4_0 = 1, x3_0*x4_5 = 1, b*x1_6*x4_2 = 1, -5*sgm*x2_4*x4_1 = 1, -35*sgm*x2_4*x4_3 = 1, b*x1_7*x4_0 = 1, b*x1_2*x4_6 = 1, b*x1_3*x4_4 = 1, x1_3*x4_5 = 1, -alpha*x1_2 = 1, -35*sgm*x2_3*x4_4 = 1, -x2_3 = 1, x1_6*x4_0 = 1, delta*sgm*x3_1*x4_4 = 1, -x1_4 = 1, x1_1*c = 1, x1_1*x4_6 = 1, delta*sgm*x3_0*x4_3 = 1, beta*x2_3 = 1, x1_1*x4_1 = 1, -sgm*x2_0*x4_1 = 1, b*x1_3*x4_2 = 1, delta*sgm*x3_0*x4_2 = 1, x1_2*c = 1, x3_7 = 1, delta*sgm*x3_6*x4_0 = 1, b*c*x1_3 = 1, x1_8*c = 1, -x1_2 = 1, b*x1_0*x4_3 = 1, delta*sgm*x3_4*x4_3 = 1, -15*sgm*x2_2*x4_4 = 1, x3_1*x4_7 = 1, delta*sgm*x3_0*x4_7 = 1, x1_1*x4_7 = 1, x3_0*x4_7 = 1, x1_7*x4_1 = 1, x1_2*x4_3 = 1, x1_4*x4_3 = 1, x1_4*x4_5 = 1, b*x1_3*x4_0 = 1, beta*x2_5 = 1, -21*sgm*x2_2*x4_5 = 1, delta*sgm*x3_3*x4_2 = 1, -sgm*x2_0*x4_6 = 1, x1_2*x4_4 = 1, x1_1*x4_3 = 1, delta*x3_1 = 1, delta*sgm*x3_3*x4_3 = 1, b*x1_0*x4_6 = 1, x3_3*x4_2 = 1, x1_3*x4_4 = 1, -sgm*x2_0*x4_3 = 1, delta*sgm*x3_1*x4_3 = 1, x1_1*x4_4 = 1, -x1_0 = 1, -20*sgm*x2_3*x4_3 = 1, delta*sgm*x3_2*x4_3 = 1, x1_2*x4_2 = 1, -7*sgm*x4_6 = 1, x2_3 = 1, x1_9*x4_0 = 1, b*x1_7*x4_1 = 1, x3_6 = 1, b*x1_4*x4_3 = 1, -x1_6 = 1, b*x1_3*x4_5 = 1, delta*sgm*x3_4*x4_0 = 1, x1_2*x4_0 = 1, b*x1_0*x4_1 = 1, b*x1_1*x4_0 = 1, b*x1_2*x4_2 = 1, b*c*x1_2 = 1, x3_2*x4_6 = 1, x3_2*x4_3 = 1, delta*sgm*x3_1*x4_1 = 1, beta*x2_4 = 1, delta*sgm*x3_7*x4_0 = 1, b*c*x1_4 = 1, b*x1_2*x4_0 = 1, -x1_8 = 1, x3_6*x4_1 = 1, delta*sgm*x3_5*x4_2 = 1, x3_1*x4_3 = 1, delta*sgm*x3_1*x4_2 = 1, x2_4 = 1, -15*sgm*x2_4*x4_2 = 1, -x2_4 = 1, delta*x3_5 = 1, x1_5*x4_4 = 1, -7*sgm*x2_6*x4_1 = 1, x1_5*x4_3 = 1, -alpha*x1_5 = 1, delta*sgm*x3_5*x4_1 = 1, x3_3*x4_5 = 1, -sgm*x2_0*x4_5 = 1, b*c*x1_5 = 1, -6*sgm*x2_5*x4_1 = 1, -x1_7 = 1, x1_5*x4_2 = 1, x3_1*x4_6 = 1, b*x1_0*x4_7 = 1, b*x1_1*x4_2 = 1, x3_0*x4_6 = 1, x1_3*c = 1, x1_4*x4_2 = 1, x1_5*x4_0 = 1, b*x1_4*x4_2 = 1, x3_2*x4_2 = 1, x3_1*x4_2 = 1, x3_4*x4_2 = 1, x3_4*x4_1 = 1, x1_1*x4_2 = 1, x1_6*x4_3 = 1, x1_9*c = 1, x3_3 = 1, beta*x2_0 = 1, x3_0*x4_8 = 1, b*x1_3*x4_1 = 1, -alpha*x1_0 = 1, x1_8*x4_0 = 1, x3_0*x4_1 = 1, z_aux*x3_0*c = 1, beta = 1, -sgm*x2_0*x4_4 = 1, b*x1_5*x4_1 = 1, b*x1_1*x4_7 = 1, -sgm*x2_0*x4_2 = 1, -x4_0*sgm = 1, b*x1_4*x4_4 = 1, x1_2*x4_6 = 1, b*x1_5*x4_2 = 1, x1_3*x4_3 = 1, b*x1_0*x4_8 = 1, delta*sgm*x3_3*x4_1 = 1, x1_3*x4_2 = 1, -x2_0 = 1, x1_3*x4_0 = 1, -sgm*x2_7*x4_0 = 1, b*x1_1*x4_6 = 1, delta*sgm*x3_0*x4_4 = 1, b*x1_1*x4_1 = 1, x3_5*x4_1 = 1, b*x1_1*x4_4 = 1, x3_1*x4_1 = 1, x3_3*x4_1 = 1, -sgm*x2_6*x4_0 = 1, delta*x3_4 = 1, x1_2*x4_1 = 1, delta*sgm*x3_0*x4_6 = 1, x1_5*c = 1, x1_6*x4_2 = 1, b*x1_1*x4_5 = 1, delta*sgm*x3_1*x4_6 = 1, b*c*x1_1 = 1, x3_5 = 1, b*c*x1_6 = 1, x3_2*x4_5 = 1, x3_2*x4_4 = 1, delta*sgm*x3_1*x4_0 = 1, x1_8*x4_1 = 1, -10*sgm*x2_2*x4_3 = 1, delta*sgm*x3_4*x4_2 = 1, -alpha*x1_4 = 1, delta*x3_3 = 1, x1_4*x4_1 = 1, x1_4*c = 1, delta*sgm*x3_2*x4_5 = 1, x1_1*x4_8 = 1, -2297306029824391595558031449960236824838226366293027466419041173826341843762215900541735527866801767603403013545732635182453132747756655164101319318982108026958082169588469820070381493073453267729147075450325641089459981159602020791921907946395140784584785846618023368797774669361300450016085425238991483784644280864603264397422678950382311065320364735644717259791515081650296884698530523574109896526731969865352999328482087931660360244433160265535043429420781305291702276797676789855187147807992622446327750748035447408/25057096371042314787586416995394769810497531634114477285846862462606603393646412137465103463968037118599488732942595890093010277896028612630792511439350000488397211447982301841892582503247439994451130958922556404167066405809035641254183621142085584562898685671827461576409918515625 = 1, x3_2*x4_1 = 1, -sgm*x2_4*x4_0 = 1, b*x1_4*x4_1 = 1, -x2_2 = 1, b*x1_2*x4_1 = 1, -alpha*x1_3 = 1, x3_7*x4_1 = 1, b*x1_0*x4_4 = 1, x3_4*x4_4 = 1, delta*sgm*x3_2*x4_1 = 1, delta*sgm*x3_6*x4_1 = 1, b*c*x1_0 = 1, b*x1_3*x4_3 = 1, -x1_3 = 1, -alpha*x1_6 = 1, -sgm*x2_3*x4_0 = 1, delta*sgm*x3_2*x4_0 = 1, -10*sgm*x2_3*x4_2 = 1, b*x1_4*x4_0 = 1, x1_4*x4_4 = 1, -4*sgm*x2_3*x4_1 = 1, delta*sgm*x3_0*x4_5 = 1, x1_4*x4_0 = 1, b*x1_8*x4_0 = 1, beta*x2_2 = 1, -alpha*x1_1 = 1, b*x1_0*x4_0 = 1, x1_6*x4_1 = 1, b*x1_6*x4_1 = 1, b*x1_5*x4_3 = 1, -628698479548236224015748601630649329907573300944618305839/16046202655528427505 = 1, b*c*x1_8 = 1, delta*sgm*x3_2*x4_2 = 1, -x2_0*x4_0*sgm = 1, x3_3*x4_4 = 1, b*x1_2*x4_5 = 1, -198426303481760771902847908118347140290786563780517916406904359248859766770598661528540187913895309318170427123345915897671959050144191463885043671133466030542288187273077216640906596212743105633232030201478398934824796560247623324096266415737986735940671145470163004991855033895730908647594846954957378739129216016464753863210182792951961640463129054406605340813613371143236908506279662955974472700784746523157528291815562761151435744700452128/67141237653384923587123056390913328422680243394455544664430931317767623284230837163376721676944495007067434272873916581429167276376476037769462872117056210207404360873122534886844340631068032316439885788183283524349276571290796459682408859375 = 1, x3_3*x4_3 = 1, x1_6*c = 1, -x1_1 = 1, delta*sgm*x3_5*x4_0 = 1, -4*sgm*x4_3 = 1, x3_0*x4_4 = 1, -39807190786230863282904089587560109583540687071500568010313300634951599389674272256774589331678477686786623900529631896602061279563381105548082895728572324021930414961462/27591166081327090329279230743899663049639823818727880781844072092272973197930329551816798532625 = 1, -x2_6 = 1, -2*sgm*x4_1 = 1, b*x1_0*x4_2 = 1, x1_5*x4_1 = 1, -x2_5 = 1, x3_5*x4_3 = 1, x1_7*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, -sgm*x2_2*x4_0 = 1, -sgm*x2_0*x4_7 = 1, x1_7*x4_0 = 1, b*x1_2*x4_4 = 1, x3_0*x4_2 = 1, x3_6*x4_2 = 1, -6*sgm*x4_5 = 1, x3_2 = 1]
# 282, -5.581916971
# 2
# [x2_1 = 10, gama = 43]
# [x2_1 = [4, 1, 4, 2, 2, 4, 4, 4, 4, 4], gama = [4, 4, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2]]