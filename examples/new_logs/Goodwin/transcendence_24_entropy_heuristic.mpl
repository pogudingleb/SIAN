infolevel[Groebner]:=10;
Et_hat := [5374729019900622199-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, 7176556909561991776-gama_0, gama_1, -x1_1-83746583531248777069172754215113009442291730878868528767/6000157425621993424, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama_0*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+35867032377819800592205749852561828759973801016344757868453128638868826461277559929014792441088692252567545797/989562229373151483209561983408736905999869211590272224, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_1-gama_0*x2_1-gama_1*x2_0)*x4_0-x4_1*(-delta*x3_0+gama_0*x2_0))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama_0*x2_0+x3_1, -x1_3-7680576074554472301944335841688413261591522706943055055807898510954114353893757200373113408787722965372360052173804391403417144342990231587135669179855600442091297/81600642811505194936814869042177135233662369791410993602857744010225131441362472014180512, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0)*x4_0+(2*delta*x3_1-2*gama_0*x2_1-2*gama_1*x2_0)*x4_1+x4_2*(delta*x3_0-gama_0*x2_0))*sgm+2*x3_1*x4_2+x3_0*x4_3+x3_2*x4_1, gama_2, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama_0*x2_1-gama_1*x2_0+x3_2, -x1_4+3316217107003235745075779042551833562995278166969601169079518366433215678190300217198501938212643609047637581217274932000206590519977826153388259032201312882024645850005662634136908092325453597339822195455282501381311/13457799236069984583175020189255640267063142958749812985740134871217602027941961936112863264318348045121309313226310832254912, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0)*x4_0+(3*delta*x3_2-3*gama_0*x2_2-6*gama_1*x2_1-3*gama_2*x2_0)*x4_1+(-x2_0*x4_3-3*x2_1*x4_2)*gama_0+(x3_0*x4_3+3*x3_1*x4_2)*delta-3*x4_2*gama_1*x2_0)*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, gama_3, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0+x3_3, -x1_5-1052916612609816343999406002582964523125420370879619190782776296017176575181403319885389182796104535953210396006654276512757634645279781152241353504412942058482999894852411037539752871629436236689964147447798984782744166032632838435650121977693782631450442408082850238678893973082912989/69359272962748538406881440091319692730388291827895057929867694921346435792746121766684561194419845797818737305337993828819692427793560379431704040451607940916, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0)*x4_0+(4*delta*x3_3-4*gama_0*x2_3-12*gama_1*x2_2-12*gama_2*x2_1-4*gama_3*x2_0)*x4_1+(-x2_0*x4_4-4*x2_1*x4_3-6*x2_2*x4_2)*gama_0+(x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2)*delta+(-4*gama_1*x4_3-6*gama_2*x4_2)*x2_0-12*x4_2*gama_1*x2_1)*sgm+4*x3_3*x4_2+x3_4*x4_1+6*x3_2*x4_3+4*x3_1*x4_4+x3_0*x4_5, gama_4, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0+x3_4, -x1_6+8275349353324654015384530015960640188223161296846712460188421440407315193957044099652999342875284138513160596618171678550126325493076366697380778826523868957667000637349261893624900170958626679830558272489276838876166276800555580775303766872210037254645801364767088854345327851472250268132497681255142934577706888883784911201787964454607074507640777834331145/91511354668969820383822117902272598818823036281183756950096734219570577355225341171088929653736701759169348436225296614053302927189504867605604802915580687397933238187337574503786246778407060928, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0)*x4_0+(5*delta*x3_4-5*gama_0*x2_4-20*gama_1*x2_3-30*gama_2*x2_2-20*gama_3*x2_1-5*gama_4*x2_0)*x4_1+(-x2_0*x4_5-5*x2_1*x4_4-10*x2_2*x4_3-10*x2_3*x4_2)*gama_0+(x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2)*delta+(-5*gama_1*x4_4-10*gama_2*x4_3-10*gama_3*x4_2)*x2_0+(-30*gama_1*x2_2-30*gama_2*x2_1)*x4_2-20*x4_3*gama_1*x2_1)*sgm+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+5*x3_1*x4_5+x3_0*x4_6, gama_5, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0+x3_5, -x1_7-451927101457289529451443156517764773281717152847753324254039863580646033719381441636574293011499035308852404541816094714683293139748536930052497601271025742394437241612575961088182901737649270422910816182061323585609425120530599694565183098890234933127233412978322393829903464743662515549745739182941139608348120565511970386551285391732143229337349013844564452580153292467242329752702438602326549160034931834036429380558068100289/7546150352029773496070513995277342457687700272470641038834455384122822452399120830855795685441372972821023565400077323703865818588751051603382844583211821508742848909358795597310691917676824007701905560397250567095706126170707264, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0)*x4_0+(6*delta*x3_5-6*gama_0*x2_5-30*gama_1*x2_4-60*gama_2*x2_3-60*gama_3*x2_2-30*gama_4*x2_1-6*gama_5*x2_0)*x4_1+(-x2_0*x4_6-6*x2_1*x4_5-15*x2_2*x4_4-20*x2_3*x4_3-15*x2_4*x4_2)*gama_0+(x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2)*delta+(-6*gama_1*x4_5-15*gama_2*x4_4-20*gama_3*x4_3-15*gama_4*x4_2)*x2_0+(-60*gama_1*x2_3-90*gama_2*x2_2-60*gama_3*x2_1)*x4_2+(-30*x2_1*x4_4-60*x2_2*x4_3)*gama_1-60*gama_2*x2_1*x4_3)*sgm+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1+x3_0*x4_7, gama_6, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0+x3_6, -x1_8-1960460638571952898307933464725576136776039540739103850594065385076699023696045338738694819523514024753385523323609601463424743651252814989068887536975561782028966795330227044131390656013563606781052468571759756864552674302549905459837652629093232932840142660069464534403796590898601908228904440063325881113732404439354239187343432962645418338668234732598903664696831859960935337220801362139087407317772198528846485588917613950638529870451044945858818784140664017454688041849674434905651199152284891839/177790224866248502197564398785275080094727550235547291382717344208469383234279256185888543322731527479923432583473525427033620475104536555244135879378835648086823564638618842522348733075329149238177911715584954060795465798814593501440111767259277736567399611649152, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((delta*x3_7-gama_0*x2_7-7*gama_1*x2_6-21*gama_2*x2_5-35*gama_3*x2_4-35*gama_4*x2_3-21*gama_5*x2_2-7*gama_6*x2_1-gama_7*x2_0)*x4_0+(7*delta*x3_6-7*gama_0*x2_6-42*gama_1*x2_5-105*gama_2*x2_4-140*gama_3*x2_3-105*gama_4*x2_2-42*gama_5*x2_1-7*gama_6*x2_0)*x4_1+(-x2_0*x4_7-7*x2_1*x4_6-21*x2_2*x4_5-35*x2_3*x4_4-35*x2_4*x4_3-21*x2_5*x4_2)*gama_0+(x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2)*delta+(-7*gama_1*x4_6-21*gama_2*x4_5-35*gama_3*x4_4-35*gama_4*x4_3-21*gama_5*x4_2)*x2_0+(-105*gama_1*x2_4-210*gama_2*x2_3-210*gama_3*x2_2-105*gama_4*x2_1)*x4_2+(-42*x2_1*x4_5-105*x2_2*x4_4-140*x2_3*x4_3)*gama_1+(-105*gama_2*x4_4-140*gama_3*x4_3)*x2_1-210*gama_2*x2_2*x4_3)*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, gama_7, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0+x3_7, -x1_9+942585438615728894653203517210470384026531514154792620178337055057606246636042955282293568288073073133417504120011873960321634890916896394885940729709444609259549633544281739315313742227481335495378560656140937060089032933174969957143186234511659522100757791650249562832143708835751162310069191911208496645463491314606608785070482382010001576587120182376828072534039966397403858147698190163527519644480946937632354217131683023653298172048277101648865368287934484173374493869116487863803878143256899806055828227417927896479480776750917120928074673782565288827807292527624131/3665205735438702049786963522621772962389917069107255680923036827275581082770788294247060302295914355201835217574258575732199929607408808104134478148466655324013652266751216390550656218594048661694597031937425394500233440603171050035567023045443909174108212947244573075910230265586097201361047598944, -gama_1, -gama_2, -gama_3, -gama_4, -gama_5, -gama_6, -gama_7, -8705455976869076300484093872274584907000234230651039116803488708680784041413046240793527269614955217911769379002144042673164345606116592474929789200336098652778784026333346817249223699273906804395406532946663314112867106796669410229283947798778628007106584373390511971007511816442005259880278903989994333/277437091850994153627525760365278770921553167311580231719470779685385743170984487066738244777679383191274949221351975315278769711174241517726816161806431763664-x2_6, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, gama_7, x4_7, x3_7, x2_7, x1_7, gama_6, x4_6, x3_6, x2_6, x1_6, gama_5, x4_5, x3_5, x2_5, x1_5, gama_4, x4_4, x3_4, x2_4, x1_4, gama_3, x4_3, x3_3, x2_3, x1_3, gama_2, x4_2, x3_2, x2_2, x1_2, gama_1, x4_1, x3_1, x2_1, x1_1, gama_0, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [gama, x2_6] 161
# 274, [1 = 7, -1 = 3, x1_2*x4_0 = 1, delta*x3_6 = 1, x1_2*x4_6 = 1, b*x1_3*x4_4 = 1, -x1_2 = 1, x3_6*x4_1 = 1, -3*sgm*x2_2*x4_1 = 1, b*x1_1*x4_2 = 1, -7*sgm*x2_1*x4_6 = 1, x3_1*x4_4 = 1, -sgm*x2_7*x4_0 = 1, x3_1*x4_7 = 1, b*c*x1_2 = 1, x1_5*x4_4 = 1, b*x1_5*x4_2 = 1, -5*sgm*x2_4*x4_1 = 1, x3_5*x4_2 = 1, x1_2*x4_4 = 1, delta*sgm*x3_1*x4_2 = 1, x1_2*x4_1 = 1, x1_1*x4_1 = 1, beta*x2_5 = 1, -1770977435351686016068346260438129682862046220309700163892609663766343204319395223811939121160074522582872295441873426109538757450368379806404807167743148460317914365946818078538984915278334890011404201915959662864877602272636609902307352381568554453912276230703458199280725004204607266229378349139204218727453582592021384244095365673331361848032317200583325478985536858394469179685595507363555448962462921753147816407357690229091755445755505/10293341692123280953159287302430384864146382227922150074030101361461568306292542099178095888330215093823260791239484669835274910692914088715490119593213342335562262684153258532947134992888627884734744720024482578963723021800411331065552 = 1, -x1_7 = 1, -sgm*x2_0*x4_4 = 1, x1_5*x4_2 = 1, -5*sgm*x2_1*x4_4 = 1, x1_8*x4_0 = 1, x1_1*x4_4 = 1, -sgm*x2_4*x4_0 = 1, delta*sgm*x3_3*x4_4 = 1, x1_1*x4_3 = 1, b*x1_4*x4_2 = 1, b*x1_1*x4_3 = 1, -15*sgm*x2_2*x4_4 = 1, delta*sgm*x3_0*x4_7 = 1, x3_0*x4_5 = 1, b*x1_4*x4_1 = 1, delta*sgm*x3_0*x4_1 = 1, x1_1*x4_5 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_3*x4_2 = 1, delta*sgm*x3_0*x4_0 = 1, -alpha*x1_6 = 1, delta*sgm*x3_2*x4_4 = 1, x1_6*x4_3 = 1, delta*sgm*x3_1*x4_1 = 1, -7*sgm*x4_1 = 1, x1_9*x4_0 = 1, b*x1_0*x4_6 = 1, x1_7*c = 1, -sgm*x2_0*x4_6 = 1, x3_2*x4_3 = 1, -sgm*x2_0*x4_7 = 1, -sgm*x2_5*x4_0 = 1, delta*sgm*x3_2*x4_3 = 1, -alpha*x1_3 = 1, x3_0*x4_3 = 1, x1_1*c = 1, b*x1_1*x4_0 = 1, delta*sgm*x3_7*x4_0 = 1, x3_4*x4_3 = 1, x1_3*x4_3 = 1, x1_7*x4_2 = 1, x3_1*x4_1 = 1, delta*sgm*x3_4*x4_0 = 1, -21*sgm*x2_2*x4_5 = 1, b*x1_2*x4_6 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_4 = 1, -10*sgm*x2_3*x4_2 = 1, b*x1_0*x4_2 = 1, b*c*x1_8 = 1, delta*x3_4 = 1, x1_7*x4_0 = 1, -x2_1 = 1, x1_3*x4_2 = 1, x1_2*x4_5 = 1, delta*sgm*x3_0*x4_3 = 1, x1_4*c = 1, b*x1_0*x4_3 = 1, b*x1_0*x4_8 = 1, x1_3*x4_4 = 1, delta*x3_1 = 1, delta*sgm*x3_2*x4_0 = 1, b*x1_1*x4_6 = 1, -alpha*x1_4 = 1, x3_2 = 1, -x1_5 = 1, x1_3*x4_0 = 1, x3_0*x4_8 = 1, x1_8*x4_1 = 1, b*x1_5*x4_3 = 1, beta*x2_1 = 1, -21*sgm*x2_5*x4_2 = 1, x3_3*x4_4 = 1, b*x1_0*x4_1 = 1, beta*x2_4 = 1, delta*sgm*x3_4*x4_3 = 1, x1_3*x4_5 = 1, -272445264882798521927162215417260549350692976219005930595/7501903452127771988 = 1, x3_4 = 1, -sgm*x2_3*x4_0 = 1, b*x1_1*x4_5 = 1, -3*sgm*x2_1*x4_2 = 1, x1_4*x4_3 = 1, b*x1_5*x4_0 = 1, x1_8*c = 1, x1_2*x4_3 = 1, delta*sgm*x3_6*x4_0 = 1, x1_2*c = 1, x1_1*x4_6 = 1, b*x1_2*x4_3 = 1, b*x1_3*x4_3 = 1, x1_1*x4_0 = 1, x3_3 = 1, delta*sgm*x3_0*x4_2 = 1, delta*sgm*x3_3*x4_3 = 1, b*x1_2*x4_2 = 1, x3_2*x4_5 = 1, b*c*x1_7 = 1, x1_5*x4_1 = 1, b*x1_6*x4_2 = 1, -20*sgm*x2_3*x4_3 = 1, z_aux*x3_0*x4_0 = 1, delta*sgm*x3_1*x4_3 = 1, b*x1_6*x4_0 = 1, -x2_0*x4_0*sgm = 1, x1_5*x4_3 = 1, b*x1_0*x4_0 = 1, delta*sgm*x3_5*x4_2 = 1, b*c*x1_1 = 1, delta*sgm*x3_1*x4_6 = 1, x1_5*c = 1, x1_3*x4_1 = 1, -x2_2 = 1, -sgm*x2_0*x4_3 = 1, x1_2*x4_7 = 1, b*x1_7*x4_0 = 1, x3_6 = 1, -4*sgm*x2_3*x4_1 = 1, beta*x2_0 = 1, b*x1_4*x4_3 = 1, delta*x3_2 = 1, -35*sgm*x2_3*x4_4 = 1, x1_3*c = 1, -sgm*x2_1*x4_0 = 1, x3_2*x4_1 = 1, b*x1_1*x4_1 = 1, x3_5*x4_3 = 1, delta*sgm*x3_3*x4_0 = 1, -sgm*x2_0*x4_2 = 1, x1_1*x4_8 = 1, -x2_5 = 1, x2_4 = 1, -35*sgm*x2_4*x4_3 = 1, delta*sgm*x3_1*x4_0 = 1, x1_1*x4_2 = 1, -x1_8 = 1, x3_4*x4_4 = 1, -6*sgm*x2_1*x4_5 = 1, x3_2*x4_4 = 1, x1_4*x4_0 = 1, x3_5 = 1, b*x1_0*x4_7 = 1, delta*x3_3 = 1, -x1_1 = 1, -x2_0 = 1, -x1_6 = 1, delta*sgm*x3_0*x4_5 = 1, b*x1_3*x4_0 = 1, -x2_3 = 1, x3_1*x4_6 = 1, delta*x3_0 = 1, b*x1_2*x4_5 = 1, x1_5*x4_0 = 1, x1_3*x4_6 = 1, -sgm*x2_0*x4_1 = 1, -x4_0*sgm = 1, x3_0*x4_4 = 1, x3_1*x4_5 = 1, x3_3*x4_5 = 1, x3_6*x4_2 = 1, b*x1_1*x4_4 = 1, delta*sgm*x3_3*x4_2 = 1, -x1_0 = 1, delta*sgm*x3_4*x4_2 = 1, -x1_4 = 1, x2_1 = 1, beta*x2_2 = 1, b*x1_4*x4_0 = 1, x3_1*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, x3_3*x4_1 = 1, b*x1_2*x4_4 = 1, -6*sgm*x2_5*x4_1 = 1, x1_4*x4_2 = 1, -alpha*x1_1 = 1, -2*sgm*x2_1*x4_1 = 1, b*x1_7*x4_1 = 1, b*c*x1_4 = 1, -4*sgm*x2_1*x4_3 = 1, x2_5 = 1, x3_4*x4_2 = 1, b*x1_5*x4_1 = 1, -6*sgm*x2_2*x4_2 = 1, x3_0*x4_7 = 1, b*x1_8*x4_0 = 1, delta*sgm*x3_6*x4_1 = 1, b*c*x1_5 = 1, x3_2*x4_2 = 1, b*x1_1*x4_7 = 1, x1_1*x4_7 = 1, x1_6*c = 1, delta*sgm*x3_2*x4_5 = 1, b*c*x1_6 = 1, b*x1_3*x4_1 = 1, -x1_3 = 1, -x2_4 = 1, delta*sgm*x3_2*x4_2 = 1, b*x1_6*x4_1 = 1, x1_6*x4_2 = 1, x3_1 = 1, x1_4*x4_5 = 1, x3_4*x4_1 = 1, -10*sgm*x2_2*x4_3 = 1, x1_6*x4_1 = 1, b*x1_2*x4_0 = 1, b*x1_4*x4_4 = 1, delta*sgm*x3_2*x4_1 = 1, beta = 1, -15*sgm*x2_4*x4_2 = 1, -alpha*x1_5 = 1, x3_0*x4_6 = 1, x1_4*x4_1 = 1, b*x1_2*x4_1 = 1, delta*x3_5 = 1, x3_2*x4_6 = 1, -sgm*x2_0*x4_5 = 1, delta*sgm*x3_3*x4_1 = 1, x3_0*x4_2 = 1, delta*sgm*x3_4*x4_1 = 1, x3_5*x4_1 = 1, x1_6*x4_0 = 1, x2_3 = 1, -alpha*x1_2 = 1, x3_7 = 1, x3_7*x4_1 = 1, x3_1*x4_3 = 1, x3_3*x4_3 = 1, delta*sgm*x3_1*x4_4 = 1, x1_2*x4_2 = 1, x2_2 = 1, -3336032149181473952864020087713184560238013976665987379171976632168321707062891453413701253431639255736689723519075308913578088832994438315617585162072130475330970534554109591586746796326537872800649354238067236193744601802894319479916635881922290489993337128746607006849023334969624596406736855010544450001221912630850208156451707132759842257636625339198715799064515/3075687755315850387175662323222007972494884287560494128503145195310479376820213705931964526464238216382270986848915206535415505804633662924467593798277116006413740347229059296594177936399730660633988 = 1, b*c*x1_3 = 1, -4556864563946229623878477120781517946238780291189388928946279646076171458135587639108897691164995222174383671048060933787619955927621795191956119414764057246930420885/2625728294855175635851655717123462957209652192301153216647477775470267752137313425729678394 = 1, b*x1_3*x4_5 = 1, z_aux*x3_0*c = 1, beta*x2_3 = 1, b*c*x1_0 = 1, x3_0*x4_1 = 1, x1_7*x4_1 = 1, b*x1_0*x4_5 = 1, -alpha*x1_0 = 1, x2_7 = 1, x3_3*x4_2 = 1, delta*sgm*x3_5*x4_1 = 1, -x1_9 = 1, delta*sgm*x3_0*x4_6 = 1, -sgm*x2_2*x4_0 = 1, b*x1_0*x4_4 = 1, x1_9*c = 1]
# 282, -5.581916971
# 2
# [x2_6 = 2, gama = 43]
# [x2_6 = [1, 2], gama = [4, 4, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2]]