infolevel[Groebner]:=10;
Et_hat := [4648889956923874854-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-71664364404343964089859083274545550927750850451073624159/2162973378189681420, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+812161456021236309074648696871911226318078899322518099803559142562388809176876667375986584086011981866404484971/3439437977725289068461802652301526620667201836794806800, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-9204103547544367310485049639015358886349254483000738020279451122988516669709031872295982031249928336605941233502343598487078545722036407873484286470627564069850223303/5469199816282537932767283162072011819030003379235566473073537494819944563160867168969272000, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+727790412053700245871962783791408096304808658485747560014155157134999569935587463701816642827669037187422665679061082030954449468941781213456024744242338641879413333584484807629162581939246819417050500566207792524640459/60394530104724326945994967623252709123781101274096846584782776954759987178980114029633410954221697410645498862135729085770000, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5-537848802105503624265767394853625067304203432108625645221760271756090881431534283855088885283211228511653838523576756061202796856873872861673784041917080737283312807543455400193209751227501377081439253947902607327560906783407536222971188005427796175855462278773803687124094832693635807591301/2304863213762726850055406548262879321069160212797796113544850436736261870106974840641686461204052918991386202110224223647767125181862461772999713880816444539200000, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6-909849990902052300657245447230319299893221985005821952364947544186132411072675866702189095802651232461889265599097447199162588028153301223921323935092664951358891589523561764413269089476319808591500741443825240670097815336282563395676498251933218282767288868222199448100367605313694766593091277356138781119839896135277870104090195434468736162679865645632435965487/146602526888470829833859591701480822723158153397043478628148617044538551179579924896601301126542303398281191330623869694115277842022265462297662909394277755698035922721750394952148910570566782720000, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7+22083694172917052068843532567786479833239604182002294387235864286642628736239408122942973910104045658010679803248642442341789884752203089785321398522360734802162592044923590590931062863497314407070373405858678523960971928485478655360444247707458864462612603233920296067516376063825549873729870512032318843690875456836712080089901426138602408356243252497958213603588514436672726279073877506984082166043204465575679845720025666018079220383/388531749623843554427484341490420321701351860747036301521891145826448495083398945520634387228780387996452113323142190453497583268532211663836977236602103842245706778031515436144647244527508530343119266135213700871555915759486848000000, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8+107591375064267327678618877433544176433663303185308779601301197327776079801900901227873048758295863698169712791621521231002614216352360825904595328321403100066601529399882156423496528064312350751169025282463612034189387409310272639111514008768960065847726696210665524555802897915388204525147666310778532846834050584138994331465816237741102438941572087694991612424767048768103579775526455193520690999664851165100159204591432858442708942933500996511983596278869223966393977413959115465410077539108352220205447927/77227652723513525943068996676206670673938234785508635568990986174046722809833107658986811248783109440343764264978505071618289430017231070783168992824752887566035381328608192636205630586654796811798063266370441501501760012114494969046644112828069293403667304274240000000, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9-5778675981988509444298552899293427106931849151374517307905210741388036907421857676977759647699687172259752887472308825076269987668449137635163317212189795016349148359201747195590052301119634292338135630042139891781927434525764080077850191061472465314118164279547548302888400102000689787652846986260957006533959018382612441988122318549752652404759829061348007768585266794185693622224274595318766546331665001929538232320760097057102534521544633826686759106784650825117090795669340435579858489920071466383865241258360667128312868669287514577665132691616645597296792743957832517174109133/245606094264687540638519199187760948079703119332964148421852435468202405667243207476920216602831033402759789170964046950333877228976642878901495221274530310100150112453593958325612061807853743560242144346939446596802003795511733549899179637125536727676551202379237342514817148778216598268615089139200000000, -1986374753365230967965039190485744350988152111935193531622631813751176747366699039559248387698361131995971951948142542265129771643373539196379338533732822515623697229627071587467438523/45576665135687816106394026350600098491916694826963053942279479123499538026340559741410600-x2_4, 7852757533120234071362251534385226926-x3_1, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [x2_4, x3_1] 66
# 274, [1 = 8, -1 = 2, -7*gama*sgm*x2_6*x4_1 = 1, x1_2*x4_0 = 1, delta*x3_6 = 1, x1_2*x4_6 = 1, b*x1_3*x4_4 = 1, -x1_2 = 1, x3_6*x4_1 = 1, b*x1_1*x4_2 = 1, b*c*x1_2 = 1, x1_5*x4_4 = 1, b*x1_5*x4_2 = 1, x3_5*x4_2 = 1, x1_2*x4_4 = 1, x1_2*x4_1 = 1, x1_1*x4_1 = 1, beta*x2_5 = 1, -gama*sgm*x2_0*x4_1 = 1, -1770977435351686016068346260438129682862046220309700163892609663766343204319395223811939121160074522582872295441873426109538757450368379806404807167743148460317914365946818078538984915278334890011404201915959662864877602272636609902307352381568554453912276230703458199280725004204607266229378349139204218727453582592021384244095365673331361848032317200583325478985536858394469179685595507363555448962462921753147816407357690229091755445755505/10293341692123280953159287302430384864146382227922150074030101361461568306292542099178095888330215093823260791239484669835274910692914088715490119593213342335562262684153258532947134992888627884734744720024482578963723021800411331065552 = 1, -x1_7 = 1, -35*gama*x4_3*sgm = 1, x4_6*delta*sgm = 1, x4_5*delta*sgm = 1, x1_5*x4_2 = 1, -20*gama*sgm*x2_3*x4_3 = 1, -gama*sgm*x2_5*x4_0 = 1, x1_8*x4_0 = 1, x1_1*x4_4 = 1, -gama*sgm*x2_0*x4_0 = 1, beta*x2_6 = 1, delta*sgm*x3_3*x4_4 = 1, x1_1*x4_3 = 1, b*x1_4*x4_2 = 1, -gama*sgm*x2_0*x4_2 = 1, b*x1_1*x4_3 = 1, delta*sgm*x3_0*x4_7 = 1, x3_0*x4_5 = 1, -2*gama*sgm*x2_1*x4_1 = 1, b*x1_4*x4_1 = 1, -gama*sgm*x2_1*x4_0 = 1, delta*sgm*x3_0*x4_1 = 1, x1_1*x4_5 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_3*x4_2 = 1, -gama*sgm*x2_3*x4_0 = 1, delta*sgm*x3_0*x4_0 = 1, -alpha*x1_6 = 1, delta*sgm*x3_2*x4_4 = 1, x1_6*x4_3 = 1, x1_9*x4_0 = 1, b*x1_0*x4_6 = 1, x1_7*c = 1, -21*gama*sgm*x2_5*x4_2 = 1, -15*gama*sgm*x2_2*x4_4 = 1, x3_2*x4_3 = 1, x4_3 = 1, delta*sgm*x3_2*x4_3 = 1, -alpha*x1_3 = 1, x4_2 = 1, x3_0*x4_3 = 1, x1_1*c = 1, b*x1_1*x4_0 = 1, delta*sgm*x3_7*x4_0 = 1, x3_4*x4_3 = 1, x1_3*x4_3 = 1, x1_7*x4_2 = 1, delta*x4_0*sgm = 1, delta*sgm*x3_4*x4_0 = 1, b*x1_2*x4_6 = 1, x1_4*x4_4 = 1, x4_4*delta*sgm = 1, delta*sgm*x3_0*x4_4 = 1, x4_2*delta*sgm = 1, -3*gama*sgm*x2_2*x4_1 = 1, -gama*sgm*x2_6*x4_0 = 1, b*x1_0*x4_2 = 1, b*c*x1_8 = 1, delta*x3_4 = 1, x4_3*delta*sgm = 1, -gama*sgm*x2_2*x4_0 = 1, -15*gama*x4_2*sgm = 1, x1_7*x4_0 = 1, x1_3*x4_2 = 1, x1_2*x4_5 = 1, delta*sgm*x3_0*x4_3 = 1, x1_4*c = 1, b*x1_0*x4_3 = 1, b*x1_0*x4_8 = 1, x1_3*x4_4 = 1, delta*sgm*x3_2*x4_0 = 1, b*x1_1*x4_6 = 1, -alpha*x1_4 = 1, x3_2 = 1, -x1_5 = 1, x1_3*x4_0 = 1, x3_0*x4_8 = 1, x1_8*x4_1 = 1, b*x1_5*x4_3 = 1, beta*x2_1 = 1, delta*x4_1*sgm = 1, -gama*sgm*x2_7*x4_0 = 1, x3_3*x4_4 = 1, b*x1_0*x4_1 = 1, delta*sgm*x3_4*x4_3 = 1, -gama*x2_3 = 1, x1_3*x4_5 = 1, x4_4 = 1, -3*gama*sgm*x2_1*x4_2 = 1, -272445264882798521927162215417260549350692976219005930595/7501903452127771988 = 1, x3_4 = 1, delta = 1, b*x1_1*x4_5 = 1, x1_4*x4_3 = 1, b*x1_5*x4_0 = 1, x1_8*c = 1, x1_2*x4_3 = 1, delta*sgm*x3_6*x4_0 = 1, x1_2*c = 1, x1_1*x4_6 = 1, b*x1_2*x4_3 = 1, b*x1_3*x4_3 = 1, x1_1*x4_0 = 1, x3_3 = 1, delta*sgm*x3_0*x4_2 = 1, delta*sgm*x3_3*x4_3 = 1, -6*gama*sgm*x2_5*x4_1 = 1, b*x1_2*x4_2 = 1, x3_2*x4_5 = 1, b*c*x1_7 = 1, x1_5*x4_1 = 1, b*x1_6*x4_2 = 1, z_aux*x3_0*x4_0 = 1, x4_5 = 1, b*x1_6*x4_0 = 1, x1_5*x4_3 = 1, b*x1_0*x4_0 = 1, delta*sgm*x3_5*x4_2 = 1, -7*gama*sgm*x2_1*x4_6 = 1, b*c*x1_1 = 1, x1_5*c = 1, x1_3*x4_1 = 1, x1_2*x4_7 = 1, b*x1_7*x4_0 = 1, x3_6 = 1, -5*gama*sgm*x2_1*x4_4 = 1, beta*x2_0 = 1, b*x1_4*x4_3 = 1, delta*x3_2 = 1, x1_3*c = 1, x3_2*x4_1 = 1, b*x1_1*x4_1 = 1, x3_5*x4_3 = 1, delta*sgm*x3_3*x4_0 = 1, x1_1*x4_8 = 1, x4_1 = 1, x1_1*x4_2 = 1, -x1_8 = 1, x3_4*x4_4 = 1, x3_2*x4_4 = 1, x1_4*x4_0 = 1, x3_5 = 1, b*x1_0*x4_7 = 1, -gama = 1, delta*x3_3 = 1, -gama*x2_1 = 1, -x1_1 = 1, -x1_6 = 1, -gama*sgm*x2_0*x4_7 = 1, delta*sgm*x3_0*x4_5 = 1, b*x1_3*x4_0 = 1, -gama*sgm*x2_0*x4_5 = 1, delta*x3_0 = 1, b*x1_2*x4_5 = 1, x1_5*x4_0 = 1, x1_3*x4_6 = 1, x3_0*x4_4 = 1, x3_3*x4_5 = 1, x3_6*x4_2 = 1, b*x1_1*x4_4 = 1, delta*sgm*x3_3*x4_2 = 1, -x1_0 = 1, delta*sgm*x3_4*x4_2 = 1, -gama*sgm*x2_0*x4_3 = 1, -x1_4 = 1, x2_1 = 1, beta*x2_2 = 1, b*x1_4*x4_0 = 1, x3_3*x4_1 = 1, b*x1_2*x4_4 = 1, -gama*sgm*x4_0 = 1, x1_4*x4_2 = 1, -alpha*x1_1 = 1, b*x1_7*x4_1 = 1, b*c*x1_4 = 1, -4*gama*sgm*x2_1*x4_3 = 1, x2_5 = 1, x3_4*x4_2 = 1, b*x1_5*x4_1 = 1, -gama*x2_5 = 1, x3_0*x4_7 = 1, b*x1_8*x4_0 = 1, delta*sgm*x3_6*x4_1 = 1, b*c*x1_5 = 1, x3_2*x4_2 = 1, b*x1_1*x4_7 = 1, x1_1*x4_7 = 1, x1_6*c = 1, -gama*x2_2 = 1, -gama*x2_6 = 1, delta*sgm*x3_2*x4_5 = 1, b*c*x1_6 = 1, -4*gama*sgm*x2_3*x4_1 = 1, b*x1_3*x4_1 = 1, -x1_3 = 1, -35*gama*sgm*x2_3*x4_4 = 1, delta*sgm*x3_2*x4_2 = 1, b*x1_6*x4_1 = 1, x1_6*x4_2 = 1, x1_4*x4_5 = 1, -gama*sgm*x2_0*x4_6 = 1, x3_4*x4_1 = 1, x1_6*x4_1 = 1, b*x1_2*x4_0 = 1, b*x1_4*x4_4 = 1, -6*gama*sgm*x2_2*x4_2 = 1, delta*sgm*x3_2*x4_1 = 1, beta = 1, -gama*x2_0 = 1, -alpha*x1_5 = 1, x3_0*x4_6 = 1, x1_4*x4_1 = 1, b*x1_2*x4_1 = 1, delta*x3_5 = 1, x3_2*x4_6 = 1, delta*sgm*x3_3*x4_1 = 1, x3_0*x4_2 = 1, delta*sgm*x3_4*x4_1 = 1, x3_5*x4_1 = 1, x1_6*x4_0 = 1, x2_3 = 1, -alpha*x1_2 = 1, x3_7 = 1, x3_7*x4_1 = 1, x4_6 = 1, x3_3*x4_3 = 1, x1_2*x4_2 = 1, -6*gama*sgm*x2_1*x4_5 = 1, -10*gama*sgm*x2_3*x4_2 = 1, x2_2 = 1, -3336032149181473952864020087713184560238013976665987379171976632168321707062891453413701253431639255736689723519075308913578088832994438315617585162072130475330970534554109591586746796326537872800649354238067236193744601802894319479916635881922290489993337128746607006849023334969624596406736855010544450001221912630850208156451707132759842257636625339198715799064515/3075687755315850387175662323222007972494884287560494128503145195310479376820213705931964526464238216382270986848915206535415505804633662924467593798277116006413740347229059296594177936399730660633988 = 1, b*c*x1_3 = 1, -4556864563946229623878477120781517946238780291189388928946279646076171458135587639108897691164995222174383671048060933787619955927621795191956119414764057246930420885/2625728294855175635851655717123462957209652192301153216647477775470267752137313425729678394 = 1, -gama*sgm*x2_0*x4_4 = 1, b*x1_3*x4_5 = 1, -5*gama*x4_1*sgm = 1, z_aux*x3_0*c = 1, beta*x2_3 = 1, b*c*x1_0 = 1, x3_0*x4_1 = 1, x1_7*x4_1 = 1, b*x1_0*x4_5 = 1, -alpha*x1_0 = 1, x2_6 = 1, -10*gama*sgm*x2_2*x4_3 = 1, x2_7 = 1, x3_3*x4_2 = 1, delta*sgm*x3_5*x4_1 = 1, -x1_9 = 1, x4_7 = 1, delta*sgm*x3_0*x4_6 = 1, -21*gama*sgm*x2_2*x4_5 = 1, b*x1_0*x4_4 = 1, x1_9*c = 1]
# 282, -5.577999884
# 2
# [x3_1 = 16, x2_4 = 7]
# [x3_1 = [4, 2, 1, 4, 2, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2], x2_4 = [4, 1, 4, 2, 2, 4, 4]]