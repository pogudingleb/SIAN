infolevel[Groebner]:=10;
Et_hat := [275649100227622010898-w_0, c*q*w_0*y_0-c*w_0*x_0*y_0+b*w_0+w_1, 35573783440295053488-z_0, -c*q*w_0*y_0+h*z_0+z_1, 156259519831120230897-k_0, k_1, -w_1-4096472762452723496142669516920379713896529958245727183770566539173311227143605500, ((-x_1*y_0+y_1*(q-x_0))*w_0+w_1*y_0*(q-x_0))*c+b*w_1+w_2, beta*v_0*x_0+d*x_0-lm+x_1, -beta*v_0*x_0+a*y_0+y_1, 6382559301896461221169525643653777266803802109578917302312486721166670381097919564-z_1, -q*(w_0*y_1+w_1*y_0)*c+h*z_1+z_2, -w_2+60878446835704426738270620928057148712979400016406169020873909955082032942216823121128653624492468425736833031063867262372110048104953556009344, ((q*y_2-x_0*y_2-2*x_1*y_1-x_2*y_0)*w_0+(-2*w_1*x_1+w_2*(q-x_0))*y_0+2*y_1*w_1*(q-x_0))*c+b*w_2+w_3, beta*x_0*v_1+(beta*v_0+d)*x_1+x_2, -beta*v_0*x_1-beta*v_1*x_0+a*y_1+y_2, -k_0*y_0+u*v_0+v_1, -94852405878951273679114352079266105386781139275169342106617694766428344150131496405777292547183734025445268455092532992724654946323832139253720-z_2, -q*(w_0*y_2+2*w_1*y_1+w_2*y_0)*c+h*z_2+z_3, -w_3-904725968911031569308190328142536707614176237763211961246390412284793159171561833250385201441473070866038849929207397941610588116435127123319285359088703814849793596871323655078153897300444220278881450004, ((q*y_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0)*w_0+(q*w_3-3*w_1*x_2-3*w_2*x_1-w_3*x_0)*y_0+(3*q*y_2-3*x_0*y_2-6*x_1*y_1)*w_1+3*y_1*w_2*(q-x_0))*c+b*w_3+w_4, (v_0*x_2+2*v_1*x_1+v_2*x_0)*beta+d*x_2+x_3, (-v_0*x_2-2*v_1*x_1-v_2*x_0)*beta+a*y_2+y_3, -k_0*y_1-k_1*y_0+u*v_1+v_2, 1409619319690773851281044502054104933136617849059180351453603483020334834414809633034276316957474184456180928139881781652070137864852490132957281592135168319930460282921828321641149449791626537068480080828-z_3, -3*(y_2*w_1+1/3*y_3*w_0+w_2*y_1+1/3*w_3*y_0)*q*c+h*z_3+z_4, -w_4+13445301602897465439540552408817044371368127947176887780875892869254457615634810374722081312487774867835453157617661947826235548836853771908931533735839042223662298830975354298118807922615630397235118718816119330439543644557317831432131546707501001272638810366551218, ((q*y_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0)*w_0+(q*w_4-4*w_1*x_3-6*w_2*x_2-4*w_3*x_1-w_4*x_0)*y_0+(4*q*y_3-4*x_0*y_3-12*x_1*y_2-12*x_2*y_1)*w_1+(-6*w_2*y_2-4*w_3*y_1)*x_0+(4*q*w_3-12*w_2*x_1)*y_1+6*y_2*q*w_2)*c+b*w_4+w_5, (v_0*x_3+3*v_1*x_2+3*v_2*x_1+v_3*x_0)*beta+d*x_3+x_4, (-v_0*x_3-3*v_1*x_2-3*v_2*x_1-v_3*x_0)*beta+a*y_3+y_4, -k_0*y_2-2*k_1*y_1-k_2*y_0+u*v_2+v_3, k_2, -20948615989575249960966205604275051185220962184555568392789623906903445296896045704860587027048565393040143532538168118299153492881935752029991145046337547686499924887091011395488918243748531974537272934243694628245996297063303876252851754113751220407464387126531588-z_4, -4*q*(w_3*y_1+3/2*y_2*w_2+y_3*w_1+1/4*y_4*w_0+1/4*w_4*y_0)*c+h*z_4+z_5, -w_5-199813138347810832850006890782484896377304483832256480606944661217804810916946682742415840779579064128437137627114428413131717077099763495878767179380509646936652518300552439532458392988997068858743344942214880760709414992791482155617078351626313670976996638848025072866157000723376566126827583098346525642361184230060353697944, ((q*y_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0)*w_0+(q*w_5-5*w_1*x_4-10*w_2*x_3-10*w_3*x_2-5*w_4*x_1-w_5*x_0)*y_0+(5*q*y_4-5*x_0*y_4-20*x_1*y_3-30*x_2*y_2-20*x_3*y_1)*w_1+(-10*w_2*y_3-10*w_3*y_2-5*w_4*y_1)*x_0+(5*q*w_4-30*w_2*x_2-20*w_3*x_1)*y_1+(10*w_2*y_3+10*w_3*y_2)*q-30*x_1*y_2*w_2)*c+b*w_5+w_6, (v_0*x_4+4*v_1*x_3+6*v_2*x_2+4*v_3*x_1+v_4*x_0)*beta+d*x_4+x_5, (-v_0*x_4-4*v_1*x_3-6*v_2*x_2-4*v_3*x_1-v_4*x_0)*beta+a*y_4+y_5, -k_0*y_3-3*k_1*y_2-3*k_2*y_1-k_3*y_0+u*v_3+v_4, k_3, 311321294869139927939664694021155149226395661086021605449265695105896012900439591109193551353487698538514883052773698624639173989102437386523251884863790045732960921910310849658133368900894194731942993314107427888266500629748133336024197356161568423429241563246107093675167971355554503781640920273012441063271590244314827761344-z_5, -q*(w_0*y_5+5*w_1*y_4+10*w_2*y_3+10*w_3*y_2+5*w_4*y_1+w_5*y_0)*c+h*z_5+z_6, -w_6+2969460368802547591236296609018984774876848251095140858019352157926601703541326899760568231959282179013416006338952819598458021802063015158421411559797386975185395648176843597139655026132373628166196803225894933792720353589310707584629735692985726313480728851308676333689659243481297286888041650529395009297569425449700956176008747697458548022081704638608339421889694442923277876130273462, ((q*y_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0)*w_0+(q*w_6-6*w_1*x_5-15*w_2*x_4-20*w_3*x_3-15*w_4*x_2-6*w_5*x_1-w_6*x_0)*y_0+(6*q*y_5-6*x_0*y_5-30*x_1*y_4-60*x_2*y_3-60*x_3*y_2-30*x_4*y_1)*w_1+(-15*w_2*y_4-20*w_3*y_3-15*w_4*y_2-6*w_5*y_1)*x_0+(6*q*w_5-60*w_2*x_3-60*w_3*x_2-30*w_4*x_1)*y_1+(15*w_2*y_4+20*w_3*y_3+15*w_4*y_2)*q+(-60*x_1*y_3-90*x_2*y_2)*w_2-60*x_1*y_2*w_3)*c+b*w_6+w_7, (v_0*x_5+5*v_1*x_4+10*v_2*x_3+10*v_3*x_2+5*v_4*x_1+v_5*x_0)*beta+d*x_5+x_6, (-v_0*x_5-5*v_1*x_4-10*v_2*x_3-10*v_3*x_2-5*v_4*x_1-v_5*x_0)*beta+a*y_5+y_6, -k_0*y_4-4*k_1*y_3-6*k_2*y_2-4*k_3*y_1-k_4*y_0+u*v_4+v_5, k_4, -w_7-44129705157526590927097738399102893793274335517520263391867438364292460900490251198128393035334815686703326299417166045865944128327017383406896187534957463827501755224810178433683178232364285664616562153290498651474638551029550700223234553367047422202473331397752349402092234795637418305845121919451844180459380942453427094051799446126861064334794356967286136243175406713092158524808209175662848553699070305194623592402604845565533432582990147044104, -4626603909643919202556977320345546142576606913819197313534566107495489012148346585816172142847608109555139269844459327924826613130523584913287900985272426585698510516555912270799910337415838238585855902744619123280730646729771217325171259516678557245266682712669670823180434292147274536452288848543132857088064395724273703077144476151985081940644415599980319625636580675881014866148732460-z_6, -k_1, -k_2, -k_3, -k_4, -3536076417424577992865076233351664167270170422314240990750802676269127140869142337542496892328616543155589907982291746724449099375585891316813970920175046505540466030461698453611380258178824252563494318-v_5, z_aux-1];
vars:=[w_7, z_6, y_6, x_6, w_6, z_5, y_5, x_5, w_5, v_5, z_4, y_4, x_4, w_4, v_4, k_4, z_3, y_3, x_3, w_3, v_3, k_3, z_2, y_2, x_2, w_2, v_2, k_2, z_1, y_1, x_1, w_1, v_1, k_1, z_0, y_0, x_0, w_0, v_0, k_0, z_aux, w_aux, a, b, beta, c, d, h, lm, q, u];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [k, v_5] 17
# 258, [1 = 16, v_4 = 1, a*y_0 = 1, b*w_1 = 1, -10*c*w_3*x_0*y_2 = 1, -w_0 = 1, -6*c*w_2*x_2*y_0 = 1, -2*c*w_1*x_0*y_1 = 1, -4*c*q*w_1*y_3 = 1, -15*c*w_4*x_2*y_0 = 1, -c*w_1*x_0*y_0 = 1, c*q*w_2*y_3 = 1, c*q*w_0*y_2 = 1, -5*c*w_1*x_0*y_4 = 1, beta*v_0*x_5 = 1, -c*w_0*x_0*y_4 = 1, -beta*v_0*x_3 = 1, -c*w_0*x_5*y_0 = 1, c*q*w_1*y_2 = 1, u*v_4 = 1, beta*v_3*x_2 = 1, -c*q*w_5*y_0 = 1, w_1 = 1, beta*v_3*x_1 = 1, -15*c*w_0*x_4*y_2 = 1, -3*c*w_1*x_2*y_0 = 1, -c*w_0*x_1*y_0 = 1, -4*c*w_3*x_0*y_1 = 1, -2*c*w_1*x_1*y_0 = 1, -10*c*q*w_3*y_2 = 1, c*q*w_3*y_1 = 1, h*z_2 = 1, -6*c*w_0*x_5*y_1 = 1, -3*c*w_2*x_0*y_1 = 1, -6*c*w_0*x_2*y_2 = 1, b*w_6 = 1, -12*c*w_2*x_1*y_1 = 1, y_2 = 1, -c*w_0*x_2*y_0 = 1, -30*c*w_4*x_1*y_1 = 1, c*q*w_5*y_1 = 1, -c*w_3*x_0*y_0 = 1, -3*beta*v_1*x_2 = 1, c*q*w_2*y_4 = 1, -z_5 = 1, -15*c*w_4*x_0*y_2 = 1, -3*c*w_2*x_1*y_0 = 1, b*w_2 = 1, -c*w_0*x_0*y_1 = 1, v_1 = 1, -z_3 = 1, -y_1 = 1, -w_6 = 1, -5*c*w_0*x_4*y_1 = 1, -c*q*w_0*y_4 = 1, beta*v_1*x_2 = 1, -5*c*q*w_4*y_1 = 1, -4*c*w_0*x_1*y_3 = 1, beta*v_0*x_1 = 1, -20*c*w_3*x_1*y_1 = 1, -c*w_0*x_0*y_0 = 1, -c*w_0*x_6*y_0 = 1, b*w_3 = 1, z_2 = 1, y_3 = 1, -z_0 = 1, c*q*w_5*y_0 = 1, x_6 = 1, -10*c*q*w_2*y_3 = 1, d*x_2 = 1, z_5 = 1, -c*q*w_0*y_2 = 1, -c*w_0*x_3*y_0 = 1, beta*v_2*x_2 = 1, -c*q*w_0*y_3 = 1, beta*v_1*x_1 = 1, -10*c*w_0*x_3*y_2 = 1, u*v_1 = 1, beta*v_1*x_4 = 1, -30*c*w_2*x_1*y_2 = 1, w_7 = 1, a*y_3 = 1, beta*x_0*v_1 = 1, -2*c*q*w_1*y_1 = 1, -w_2 = 1, z_4 = 1, -20*c*w_3*x_3*y_0 = 1, c*q*w_2*y_1 = 1, -c*q*w_4*y_0 = 1, beta*v_2*x_3 = 1, -6*c*w_1*x_0*y_5 = 1, -w_5 = 1, -6*c*w_2*x_0*y_2 = 1, y_1 = 1, -6*c*q*w_2*y_2 = 1, c*q*w_2*y_2 = 1, -w_7 = 1, -10*c*w_0*x_2*y_3 = 1, -3*c*w_1*x_0*y_2 = 1, -10*c*w_2*x_3*y_0 = 1, -y_4 = 1, -12*c*w_1*x_2*y_1 = 1, d*x_4 = 1, beta*v_2*x_1 = 1, -60*c*w_3*x_2*y_1 = 1, c*q*w_6*y_0 = 1, v_2 = 1, -c*q*w_0*y_0 = 1, -beta*x_0 = 1, c*q*w_1*y_5 = 1, -c*q*w_0*y_5 = 1, -c*q*w_1*y_0 = 1, a*y_1 = 1, beta*x_0 = 1, -4*c*w_1*x_3*y_0 = 1, beta*v_2*x_0 = 1, -c*q*w_0*y_1 = 1, c*q*w_0*y_1 = 1, -3*c*q*w_1*y_2 = 1, beta*v_0*x_2 = 1, -30*c*w_1*x_1*y_4 = 1, x_3 = 1, -5*c*w_0*x_1*y_4 = 1, -w_3 = 1, -5*c*q*w_1*y_4 = 1, -c*w_2*x_0*y_0 = 1, -4*c*q*w_3*y_1 = 1, -30*c*w_2*x_2*y_1 = 1, h*z_4 = 1, -c*q*w_2*y_0 = 1, -6*c*w_1*x_1*y_1 = 1, -20*c*w_1*x_3*y_1 = 1, -3*c*q*w_2*y_1 = 1, -beta*v_0*x_4 = 1, c*q*w_0*y_4 = 1, -3*c*w_0*x_1*y_2 = 1, -20*c*w_3*x_0*y_3 = 1, -5*beta*v_1*x_4 = 1, b*w_4 = 1, -z_6 = 1, c*q*w_4*y_2 = 1, -4*c*w_1*x_0*y_3 = 1, h*z_0 = 1, beta*v_0*x_3 = 1, -beta*v_0*x_1 = 1, c*q*w_0*y_6 = 1, -c*w_0*x_0*y_5 = 1, -beta*v_0*x_2 = 1, -c*w_0*x_0*y_6 = 1, c*q*w_4*y_0 = 1, -60*c*w_1*x_2*y_3 = 1, c*q*w_3*y_0 = 1, -12*c*w_1*x_1*y_2 = 1, c*q*w_0*y_3 = 1, b*w_5 = 1, -c*w_0*x_0*y_3 = 1, x_2 = 1, beta*v_1*x_3 = 1, -5*c*w_1*x_4*y_0 = 1, -c*q*w_3*y_0 = 1, z_1 = 1, -y_3 = 1, -5*c*w_4*x_1*y_0 = 1, -60*c*w_3*x_1*y_2 = 1, beta*v_0*x_4 = 1, -w_4 = 1, -w_1 = 1, c*q*w_3*y_2 = 1, -beta*x_0*v_1 = 1, beta*v_0*x_0 = 1, -30*c*w_1*x_4*y_1 = 1, x_1 = 1, beta*v_4*x_0 = 1, -10*c*w_3*x_2*y_0 = 1, -15*c*w_2*x_4*y_0 = 1, c*q*w_1*y_1 = 1, c*q*w_0*y_5 = 1, -4*beta*v_1*x_3 = 1, -3*c*w_0*x_2*y_1 = 1, -y_0 = 1, -y_2 = 1, -5*c*w_4*x_0*y_1 = 1, -z_1 = 1, -z_4 = 1, y_6 = 1, h*z_5 = 1, -10*c*w_2*x_0*y_3 = 1, -6*beta*v_2*x_2 = 1, c*q*w_0*y_0 = 1, x_4 = 1, beta*v_3*x_0 = 1, -beta*v_0*x_5 = 1, -60*c*w_1*x_3*y_2 = 1, w_5 = 1, w_2 = 1, d*x_1 = 1, -6*c*w_5*x_0*y_1 = 1, -6*c*w_1*x_5*y_0 = 1, -beta*v_3*x_0 = 1, c*q*w_1*y_3 = 1, h*z_1 = 1, h*z_3 = 1, w_3 = 1, c*q*w_3*y_3 = 1, d*x_5 = 1, v_3 = 1, -beta*v_0*x_0 = 1, -2*c*w_0*x_1*y_1 = 1, -20*c*w_0*x_3*y_3 = 1, x_5 = 1, b*w_0 = 1, -4*beta*v_3*x_1 = 1, -4*c*w_0*x_3*y_1 = 1, -15*c*w_0*x_2*y_4 = 1, beta*v_4*x_1 = 1, -c*w_4*x_0*y_0 = 1, -5*beta*v_4*x_1 = 1, y_5 = 1, a*y_2 = 1, -3*beta*v_2*x_1 = 1, -20*c*w_1*x_1*y_3 = 1, y_4 = 1, w_6 = 1, -60*c*w_2*x_3*y_1 = 1, u*v_0 = 1, -10*beta*v_2*x_3 = 1, u*v_2 = 1, -c*w_0*x_0*y_2 = 1, a*y_5 = 1, -30*c*w_1*x_2*y_2 = 1, w_4 = 1, a*y_4 = 1, c*q*w_4*y_1 = 1, -z_2 = 1, -60*c*w_2*x_1*y_3 = 1, -2*beta*v_1*x_1 = 1, -c*w_5*x_0*y_0 = 1, -4*c*w_3*x_1*y_0 = 1, -6*c*w_5*x_1*y_0 = 1, -beta*v_4*x_0 = 1, -c*w_0*x_4*y_0 = 1, d*x_0 = 1, z_3 = 1, d*x_3 = 1, u*v_3 = 1, -lm = 1, -90*c*w_2*x_2*y_2 = 1, z_6 = 1, -beta*v_2*x_0 = 1, c*q*w_1*y_0 = 1, -c*w_6*x_0*y_0 = 1, c*q*w_1*y_4 = 1, -6*c*w_0*x_1*y_5 = 1, c*q*w_2*y_0 = 1, -10*beta*v_3*x_2 = 1, z_aux = 1, -15*c*w_2*x_0*y_4 = 1, -1 = 1]
# 273, -5.446975752
# 1
# [v_5 = 3, k = 5]
# [v_5 = [3, 3, 1], k = [2, 2, 2, 2, 2]]
# [v_5 = [1, 1, 1], k = [1, 1, 1, 1, 1]]