infolevel[Groebner]:=10;
Et_hat := [216133214810640319784-w_0, c*q*w_0*y_0-c*w_0*x_0*y_0+b*w_0+w_1, 627461342895326947340-z_0, -c*q*w_0*y_0+h*z_0+z_1, 147479997310802958031-k_0, k_1, 641949061155917752926-lm_0, lm_1, -w_1+60184379324152757425395692965839971146252803156205555066210381185106171817201808, ((-x_1*y_0+y_1*(q-x_0))*w_0+w_1*y_0*(q-x_0))*c+b*w_1+w_2, beta*v_0*x_0+d*x_0-lm_0+x_1, -beta*v_0*x_0+a*y_0+y_1, 10180837867927885579754404886567911499441396448516661161321519363296442159253819300-z_1, -q*(w_0*y_1+w_1*y_0)*c+h*z_1+z_2, -w_2+16758921194999897195142533959178657853869244906495858823751343768531732981471421211164890608339100926303320727472004219566139025243682563648, ((q*y_2-x_0*y_2-2*x_1*y_1-x_2*y_0)*w_0+(-2*w_1*x_1+w_2*(q-x_0))*y_0+2*y_1*w_1*(q-x_0))*c+b*w_2+w_3, beta*x_0*v_1-lm_1+(beta*v_0+d)*x_1+x_2, -beta*v_0*x_1-beta*v_1*x_0+a*y_1+y_2, -k_0*y_0+u*v_0+v_1, 2834952548213823528079108928200388430700027604166241902441068887652874138824559544870127801170381022434318215496877337240261906973471525027500-z_2, -q*(w_0*y_2+2*w_1*y_1+w_2*y_0)*c+h*z_2+z_3, -w_3+4666683328368288057153146896692631752077385974587778780709491939109322994752244029900045481235066110358036156663493665294350039978152548811934261518576269592956098781774128576004070926291790259258880, ((q*y_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0)*w_0+(q*w_3-3*w_1*x_2-3*w_2*x_1-w_3*x_0)*y_0+(3*q*y_2-3*x_0*y_2-6*x_1*y_1)*w_1+3*y_1*w_2*(q-x_0))*c+b*w_3+w_4, (v_0*x_2+2*v_1*x_1+v_2*x_0)*beta+d*x_2-lm_2+x_3, (-v_0*x_2-2*v_1*x_1-v_2*x_0)*beta+a*y_2+y_3, lm_2, -k_0*y_1-k_1*y_0+u*v_1+v_2, 789419894009157753256255735850414104186760178106982998250703772366085318296889908133506340416910014819613637956071304602829045513042250930679845049779976783127843511018554424166342449136252158573818500-z_3, -3*(y_2*w_1+1/3*y_3*w_0+w_2*y_1+1/3*w_3*y_0)*q*c+h*z_3+z_4, -w_4+1299483005729991179123098474116660565952421450399014854651547121595595184139696900227757783101242356076212674518002210098810682136003433514704171908457041806996772597872773440383377623974004665406084148277375668391277118734771979084636728644480954689362270624, ((q*y_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0)*w_0+(q*w_4-4*w_1*x_3-6*w_2*x_2-4*w_3*x_1-w_4*x_0)*y_0+(4*q*y_3-4*x_0*y_3-12*x_1*y_2-12*x_2*y_1)*w_1+(-6*w_2*y_2-4*w_3*y_1)*x_0+(4*q*w_3-12*w_2*x_1)*y_1+6*y_2*q*w_2)*c+b*w_4+w_5, (v_0*x_3+3*v_1*x_2+3*v_2*x_1+v_3*x_0)*beta+d*x_3+x_4-lm_3, (-v_0*x_3-3*v_1*x_2-3*v_2*x_1-v_3*x_0)*beta+a*y_3+y_4, lm_3, -k_0*y_2-2*k_1*y_1-k_2*y_0+u*v_2+v_3, k_2, 219821587296080133320925597246138225149530100072648713815807445591821167513053366492578338436648381065165595372736380636544010734750693015960515571621802250370803549169258811655551573573072115768628774922320562267248751720496747083537754183420399704423913823500-z_4, -4*(w_3*y_1+3/2*y_2*w_2+y_3*w_1+1/4*y_4*w_0+1/4*w_4*y_0)*q*c+h*z_4+z_5, -w_5+361853582803848192624437569449347259957656318300138883234867251885356142577960744387281241448085189279288540386977573690680839559232118436292750541003366167143060582659257047662101969355458159849535115446142897987251879800993568394213729985567157672539763113562236361915594063005967481415423660465052317942416763255840, ((q*y_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0)*w_0+(q*w_5-5*w_1*x_4-10*w_2*x_3-10*w_3*x_2-5*w_4*x_1-w_5*x_0)*y_0+(5*q*y_4-5*x_0*y_4-20*x_1*y_3-30*x_2*y_2-20*x_3*y_1)*w_1+(-10*w_2*y_3-10*w_3*y_2-5*w_4*y_1)*x_0+(5*q*w_4-30*w_2*x_2-20*w_3*x_1)*y_1+(10*w_2*y_3+10*w_3*y_2)*q-30*x_1*y_2*w_2)*c+b*w_5+w_6, (v_0*x_4+4*v_1*x_3+6*v_2*x_2+4*v_3*x_1+v_4*x_0)*beta+d*x_4+x_5-lm_4, (-v_0*x_4-4*v_1*x_3-6*v_2*x_2-4*v_3*x_1-v_4*x_0)*beta+a*y_4+y_5, lm_4, -k_0*y_3-3*k_1*y_2-3*k_2*y_1-k_3*y_0+u*v_3+v_4, k_3, 61211442235084611495325755351920490459790291305042212372717492481885164027713747374312970349447486180700436109702226987741577710172853297076216259699232746002134664378642759247066169720528897008437847136032896711949605872630473518479910700580100109678140820959046835474877926605796049878717816125550928754365697492886500-z_5, -q*(w_0*y_5+5*w_1*y_4+10*w_2*y_3+10*w_3*y_2+5*w_4*y_1+w_5*y_0)*c+h*z_5+z_6, -w_6+100761621976292270927155003104696279763857224913335931457880002067455183226463680442684871433495653677740848209538798929954240923890814320658337861564523013631891441038618478366917478626868071216021345698423986001772197201499566794427942487417323346736476580368478320713007694180056616790622940009593767292691948160941665056895457947969451063074943322060161216096904888489474880, ((q*y_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0)*w_0+(q*w_6-6*w_1*x_5-15*w_2*x_4-20*w_3*x_3-15*w_4*x_2-6*w_5*x_1-w_6*x_0)*y_0+(6*q*y_5-6*x_0*y_5-30*x_1*y_4-60*x_2*y_3-60*x_3*y_2-30*x_4*y_1)*w_1+(-15*w_2*y_4-20*w_3*y_3-15*w_4*y_2-6*w_5*y_1)*x_0+(6*q*w_5-60*w_2*x_3-60*w_3*x_2-30*w_4*x_1)*y_1+(15*w_2*y_4+20*w_3*y_3+15*w_4*y_2)*q+(-60*x_1*y_3-90*x_2*y_2)*w_2-60*x_1*y_2*w_3)*c+b*w_6+w_7, (v_0*x_5+5*v_1*x_4+10*v_2*x_3+10*v_3*x_2+5*v_4*x_1+v_5*x_0)*beta+d*x_5+x_6-lm_5, (-v_0*x_5-5*v_1*x_4-10*v_2*x_3-10*v_3*x_2-5*v_4*x_1-v_5*x_0)*beta+a*y_5+y_6, lm_5, -k_0*y_4-4*k_1*y_3-6*k_2*y_2-4*k_3*y_1-k_4*y_0+u*v_4+v_5, k_4, -w_7+28058045977112409798139117501068441184577274479902785367021743623196840816034466295266256958704579739317796543123502274609108608010040273902483047896096445273838445319062124515023420862057013853916776771056640487508668120798460595714367008927822793005492836999760019138458853073344891251305263329086369529181075844152477684752242200056879393924842663575385766671117544786994106700202312991088615286065108240230591031627082638725977649568, 17044916773585292924091999581965986478067049636931038655497776874146310913740615104801299545856167428405807292657727125107162217714722062194576635671744299840025924884716704170487943938568189744005345790578708443269641638358797636106017119122442731873630233137973303766887349871420470517101767628993190494207855383318059272784819098782206946877161355541875062807982374619766803500-z_6, -k_1, -k_2, -k_3, -k_4, -lm_1, -lm_2, -lm_3, -lm_4, -lm_5, z_aux-1];
vars:=[w_7, z_6, y_6, x_6, w_6, lm_5, z_5, y_5, x_5, w_5, v_5, lm_4, z_4, y_4, x_4, w_4, v_4, k_4, lm_3, z_3, y_3, x_3, w_3, v_3, k_3, lm_2, z_2, y_2, x_2, w_2, v_2, k_2, lm_1, z_1, y_1, x_1, w_1, v_1, k_1, lm_0, z_0, y_0, x_0, w_0, v_0, k_0, z_aux, w_aux, a, b, beta, c, d, h, q, u];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [k, lm] 11
# 265, [1 = 8, -1 = 2, h*z_3 = 1, a*y_4 = 1, -30*c*w_1*x_1*y_4 = 1, -y_1 = 1, -3*c*w_0*x_2*y_1 = 1, -z_6 = 1, beta*v_0*x_2 = 1, h*z_5 = 1, -60*c*w_3*x_2*y_1 = 1, y_2 = 1, -6*c*w_2*x_0*y_2 = 1, -1238850523131833999692235920125556477868773187666833811120568482771192993504942636969435124262224497922864286069361616374697195553081269926755095129092181625603816676820701013546709562257982212444463663308468624912549712214904789937486399419798432795970920972608857389773268297997914605519006261200875892607715514323397348574636702830540004604507979310583346385652011268601801735497416 = 1, d*x_2 = 1, -beta*v_0*x_3 = 1, -10*c*w_3*x_0*y_2 = 1, -4*c*w_1*x_3*y_0 = 1, -w_5 = 1, -10*c*w_0*x_2*y_3 = 1, -y_3 = 1, -10*c*w_2*x_3*y_0 = 1, -5*c*q*w_1*y_4 = 1, c*q*w_5*y_0 = 1, -4*c*q*w_3*y_1 = 1, w_2 = 1, beta*v_3*x_2 = 1, c*q*w_2*y_3 = 1, -beta*x_0*v_1 = 1, c*q*w_2*y_0 = 1, h*z_0 = 1, -30*c*w_2*x_1*y_2 = 1, -beta*v_0*x_1 = 1, -c*q*w_4*y_0 = 1, -beta*v_0*x_2 = 1, -20*c*w_3*x_3*y_0 = 1, c*q*w_4*y_2 = 1, -10*c*w_3*x_2*y_0 = 1, -12*c*w_2*x_1*y_1 = 1, c*q*w_1*y_3 = 1, -w_6 = 1, -z_1 = 1, w_7 = 1, -10*c*w_0*x_3*y_2 = 1, x_6 = 1, -10*c*q*w_3*y_2 = 1, -w_7 = 1, w_5 = 1, w_3 = 1, -c*q*w_3*y_0 = 1, beta*v_0*x_0 = 1, -10*c*q*w_2*y_3 = 1, beta*v_2*x_1 = 1, a*y_1 = 1, b*w_1 = 1, c*q*w_2*y_4 = 1, beta*v_1*x_3 = 1, -30*c*w_2*x_2*y_1 = 1, -3*c*w_0*x_1*y_2 = 1, -6*c*w_1*x_0*y_5 = 1, -w_1 = 1, -20*c*w_3*x_1*y_1 = 1, -15*c*w_4*x_0*y_2 = 1, -2*beta*v_1*x_1 = 1, u*v_1 = 1, x_5 = 1, -5*c*w_4*x_0*y_1 = 1, -w_3 = 1, -10*c*w_2*x_0*y_3 = 1, -y_2 = 1, -c*q*w_0*y_0 = 1, -c*q*w_1*y_0 = 1, -6*beta*v_2*x_2 = 1, c*q*w_0*y_1 = 1, beta*v_0*x_3 = 1, -15*c*w_2*x_4*y_0 = 1, -12*c*w_1*x_1*y_2 = 1, -10591982956577674406960676710223795482278585121541118595653731851017668239845049793619418663641952992570806100098981940411406122856121024134379759707789044551475609626862887049808254061461363797250430864 = 1, -3*c*q*w_1*y_2 = 1, -y_4 = 1, beta*v_2*x_3 = 1, w_6 = 1, -20*c*w_3*x_0*y_3 = 1, z_2 = 1, y_6 = 1, -6*c*w_1*x_1*y_1 = 1, -30*c*w_1*x_4*y_1 = 1, -3*c*q*w_2*y_1 = 1, -beta*v_0*x_5 = 1, -y_0 = 1, c*q*w_1*y_4 = 1, c*q*w_1*y_1 = 1, -c*w_0*x_3*y_0 = 1, -60*c*w_1*x_2*y_3 = 1, -c*q*w_0*y_3 = 1, c*q*w_2*y_1 = 1, a*y_0 = 1, beta*v_2*x_2 = 1, y_4 = 1, -4*beta*v_1*x_3 = 1, -c*q*w_0*y_1 = 1, -c*w_0*x_0*y_6 = 1, -c*w_0*x_0*y_5 = 1, -10*beta*v_2*x_3 = 1, u*v_2 = 1, -c*q*w_0*y_5 = 1, u*v_3 = 1, -beta*v_0*x_4 = 1, -c*w_2*x_0*y_0 = 1, b*w_3 = 1, a*y_2 = 1, c*q*w_2*y_2 = 1, -5*beta*v_1*x_4 = 1, c*q*w_5*y_1 = 1, -12*c*w_1*x_2*y_1 = 1, -5*c*w_0*x_1*y_4 = 1, -160384182955497141382823800457813098826896606009047224230606208495201262522850319793793638262166142072894427665475312342224811258056808274298234867970937968055106695157950986057496934399896001764981724451943728043583837248557585624273089490857349933249291673230043661190995339790063375756735309317207621243238438142918919680 = 1, beta*v_1*x_2 = 1, c*q*w_1*y_5 = 1, -60*c*w_1*x_3*y_2 = 1, -z_0 = 1, -5*c*w_1*x_4*y_0 = 1, -3*c*w_1*x_0*y_2 = 1, -5*c*w_4*x_1*y_0 = 1, -699508523128880244923001612615625014529544595235571828887627756357009846103118048 = 1, -60*c*w_3*x_1*y_2 = 1, -c*w_0*x_2*y_0 = 1, -z_4 = 1, z_4 = 1, -c*w_3*x_0*y_0 = 1, -90*c*w_2*x_2*y_2 = 1, c*q*w_6*y_0 = 1, c*q*w_3*y_1 = 1, b*w_4 = 1, -c*w_6*x_0*y_0 = 1, beta*v_4*x_0 = 1, beta*v_2*x_0 = 1, -c*w_0*x_0*y_1 = 1, d*x_5 = 1, -c*w_5*x_0*y_0 = 1, -w_4 = 1, -2*c*q*w_1*y_1 = 1, -5403191785151578035281168156192033510041488011599520052944548122979328446614845693097117406982629088408086692928291907877174464533530506840368 = 1, -beta*v_4*x_0 = 1, -20*c*w_1*x_3*y_1 = 1, v_1 = 1, -20*c*w_0*x_3*y_3 = 1, u*v_0 = 1, d*x_1 = 1, v_4 = 1, -c*w_4*x_0*y_0 = 1, b*w_2 = 1, w_1 = 1, c*q*w_1*y_2 = 1, -c*w_0*x_1*y_0 = 1, -c*w_0*x_0*y_0 = 1, -30*c*w_1*x_2*y_2 = 1, c*q*w_1*y_0 = 1, v_3 = 1, c*q*w_3*y_2 = 1, -10*beta*v_3*x_2 = 1, -2*c*w_1*x_0*y_1 = 1, w_4 = 1, beta*v_0*x_1 = 1, c*q*w_0*y_2 = 1, z_6 = 1, -6*c*w_0*x_2*y_2 = 1, beta*v_1*x_4 = 1, c*q*w_0*y_5 = 1, -c*q*w_0*y_2 = 1, beta*v_0*x_5 = 1, -4*beta*v_3*x_1 = 1, -beta*v_3*x_0 = 1, d*x_0 = 1, v_5 = 1, -15*c*w_0*x_2*y_4 = 1, d*x_4 = 1, -5*beta*v_4*x_1 = 1, -c*w_0*x_0*y_4 = 1, -6*c*w_5*x_1*y_0 = 1, y_1 = 1, -z_3 = 1, beta*v_4*x_1 = 1, d*x_3 = 1, -6*c*w_5*x_0*y_1 = 1, u*v_4 = 1, -4*c*w_0*x_1*y_3 = 1, y_5 = 1, -c*q*w_2*y_0 = 1, -2*c*w_1*x_1*y_0 = 1, z_1 = 1, v_2 = 1, beta*v_1*x_1 = 1, beta*v_3*x_1 = 1, -20*c*w_1*x_1*y_3 = 1, z_3 = 1, -2428543007268361829797898145935981475093566844866646368399514509232827121167119206637730541505973444561807712331375742069488412264708218901026229366441883451886798621783148723794701476063569416402563930828096268625735781494111391164762975803978446612170189733726576442663412743083375205333580904038129318008688344214247980764979670478991986220856928532633541196348240594133479116710072978645613811198700053428812385068717547776784671422075400720 = 1, -6*c*w_0*x_1*y_5 = 1, b*w_0 = 1, -w_0 = 1, c*q*w_0*y_6 = 1, x_2 = 1, a*y_5 = 1, b*w_5 = 1, -60*c*w_2*x_3*y_1 = 1, -c*w_1*x_0*y_0 = 1, x_3 = 1, c*q*w_0*y_3 = 1, -beta*v_5*x_0 = 1, b*w_6 = 1, -6*c*w_2*x_2*y_0 = 1, a*y_3 = 1, -c*w_0*x_4*y_0 = 1, -3*beta*v_2*x_1 = 1, x_1 = 1, beta*v_0*x_4 = 1, h*z_2 = 1, -z_5 = 1, c*q*w_3*y_0 = 1, -c*w_0*x_6*y_0 = 1, -5*c*w_1*x_0*y_4 = 1, -c*w_0*x_0*y_3 = 1, -6*c*w_1*x_5*y_0 = 1, -3*c*w_1*x_2*y_0 = 1, -w_2 = 1, -15*c*w_4*x_2*y_0 = 1, -6*c*w_0*x_5*y_1 = 1, z_5 = 1, -c*q*w_0*y_4 = 1, -c*w_0*x_0*y_2 = 1, beta*v_5*x_0 = 1, -5*c*q*w_4*y_1 = 1, c*q*w_4*y_1 = 1, -z_2 = 1, -4*c*w_1*x_0*y_3 = 1, -c*w_0*x_5*y_0 = 1, -c*q*w_5*y_0 = 1, -beta*v_0*x_0 = 1, -2*c*w_0*x_1*y_1 = 1, c*q*w_0*y_0 = 1, -60*c*w_2*x_1*y_3 = 1, beta*x_0*v_1 = 1, -3*c*w_2*x_1*y_0 = 1, -4*c*w_0*x_3*y_1 = 1, -beta*v_2*x_0 = 1, -30*c*w_4*x_1*y_1 = 1, -81815322340113410247373836646548347083174869816821684223190817690067904156980626634559582990243776572861047930378381084997736290487009950736720452399774078602760655688831367308185001330430005011874555589058482654437908850983642811267971145276594395247095139835432 = 1, z_aux = 1, -4*c*w_3*x_0*y_1 = 1, -6*c*q*w_2*y_2 = 1, y_3 = 1, c*q*w_4*y_0 = 1, -4*c*q*w_1*y_3 = 1, beta*v_3*x_0 = 1, c*q*w_0*y_4 = 1, c*q*w_3*y_3 = 1, -4*c*w_3*x_1*y_0 = 1, -15*c*w_0*x_4*y_2 = 1, -5*c*w_0*x_4*y_1 = 1, h*z_4 = 1, -3*beta*v_1*x_2 = 1, h*z_1 = 1, x_4 = 1, -3*c*w_2*x_0*y_1 = 1, -15*c*w_2*x_0*y_4 = 1]
# 273, -5.543457778
# 2
# [lm = 1, k = 5]
# [lm = [1], k = [2, 2, 2, 2, 2]]