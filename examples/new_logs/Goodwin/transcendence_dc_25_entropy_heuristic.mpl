infolevel[Groebner]:=10;
Et_hat := [5796881316067207907-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, 7030717284584950252-gama_0, gama_1, -x1_1-168539371189470124728089775587605913195422980812122251589/14717131615698185905, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama_0*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+5690297387263435967408875616566854074294240673699137448692518650158302453363749114055188580614661746070543514833/251520265076959946505433126692628083898174483213064337835, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_1-gama_0*x2_1-gama_1*x2_0)*x4_0-x4_1*(-delta*x3_0+gama_0*x2_0))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama_0*x2_0+x3_1, -x1_3-192118222151763116789459069382348722435519459445385578626615071095794643881582004573523157027587447082964188513261198556690717007353176521113587550875764889230669409479/4298557993250847521354774160482481838763480996211300750907967736123201765499094073268125070345, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0)*x4_0+(2*delta*x3_1-2*gama_0*x2_1-2*gama_1*x2_0)*x4_1+x4_2*(delta*x3_0-gama_0*x2_0))*sgm+2*x3_1*x4_2+x3_0*x4_3+x3_2*x4_1, gama_2, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama_0*x2_1-gama_1*x2_0+x3_2, -x1_4+87178756321418121650445670467996781351438001479143876315049128037849854096092543172606154215701729029548601660227912145418077568981249801035738094709821375775444363721140828847337648055901289462020555631395298060667462657077/367318331500783713320226434401232270237174022701709565137427124724780546150212570766593163692643483401130561382195720587455657904575, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0)*x4_0+(3*delta*x3_2-3*gama_0*x2_2-6*gama_1*x2_1-3*gama_2*x2_0)*x4_1+(-x2_0*x4_3-3*x2_1*x4_2)*gama_0+(x3_0*x4_3+3*x3_1*x4_2)*delta-3*x4_2*gama_1*x2_0)*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, gama_3, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0+x3_3, -x1_5-1250166677568372695527108690183541169622876408099395892635625535985535268821764140320209508301221449933177629471963972522533220871033617738974000832749540649525049879292131997405782699507699500703425650981759618369397188071991070255883258727214599034245695021740925506617477711766045582300068930883/251103289742023135354904843596181962872014667722080311429204768981869984940535363844152309420552008566603873225423726061337243616113928814224872806712267445541268147181, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0)*x4_0+(4*delta*x3_3-4*gama_0*x2_3-12*gama_1*x2_2-12*gama_2*x2_1-4*gama_3*x2_0)*x4_1+(-x2_0*x4_4-4*x2_1*x4_3-6*x2_2*x4_2)*gama_0+(x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2)*delta+(-4*gama_1*x4_3-6*gama_2*x4_2)*x2_0-12*x4_2*gama_1*x2_1)*sgm+4*x3_3*x4_2+x3_4*x4_1+6*x3_2*x4_3+4*x3_1*x4_4+x3_0*x4_5, gama_4, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0+x3_4, -x1_6-6423364018620778245031557019973263399979606765379520449844980718666302250426672506316971376538907822162534701766954338690046542649169850287638867928760710480623440341158441169852172272371210104493049157925711442783116983149105739198509258265313205087625700529312457964005651122210142541932450265171329859008689638116714081150509164862691620474115509157798712380232587478811/107285793942079451908532008136438394718565419694784516035200034411054278177729686642663637168356476635571267461323562022121814437887269521837080238215301498158723992727723979537478726066676419662942385939175, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0)*x4_0+(5*delta*x3_4-5*gama_0*x2_4-20*gama_1*x2_3-30*gama_2*x2_2-20*gama_3*x2_1-5*gama_4*x2_0)*x4_1+(-x2_0*x4_5-5*x2_1*x4_4-10*x2_2*x4_3-10*x2_3*x4_2)*gama_0+(x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2)*delta+(-5*gama_1*x4_4-10*gama_2*x4_3-10*gama_3*x4_2)*x2_0+(-30*gama_1*x2_2-30*gama_2*x2_1)*x4_2-20*x4_3*gama_1*x2_1)*sgm+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+5*x3_1*x4_5+x3_0*x4_6, gama_5, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0+x3_5, -x1_7+68929353987298016703379953984175578592551387442151037994296780481165848478559720197473195350589739657645347153576800369034251658657232917202615036620320262058517400440944427645223960388569499523725668133791591522622901214458005116180866582785293063686408819493854505362364962692464918362326525822079704777329994237613679390837979565685931882935727940679621171959006233563311456190394501961534490597059072574207624989188045410827560704323276585327609/9167734595279615138042326025381012882402547048261127176227203503919471382286186957082823124757507732798802317193295240490794794792534438750924440083357896233557858267208654060299370223338723113639249429483047705789865794594802867593798000508625, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0)*x4_0+(6*delta*x3_5-6*gama_0*x2_5-30*gama_1*x2_4-60*gama_2*x2_3-60*gama_3*x2_2-30*gama_4*x2_1-6*gama_5*x2_0)*x4_1+(-x2_0*x4_6-6*x2_1*x4_5-15*x2_2*x4_4-20*x2_3*x4_3-15*x2_4*x4_2)*gama_0+(x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2)*delta+(-6*gama_1*x4_5-15*gama_2*x4_4-20*gama_3*x4_3-15*gama_4*x4_2)*x2_0+(-60*gama_1*x2_3-90*gama_2*x2_2-60*gama_3*x2_1)*x4_2+(-30*x2_1*x4_4-60*x2_2*x4_3)*gama_1-60*gama_2*x2_1*x4_3)*sgm+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1+x3_0*x4_7, gama_6, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0+x3_6, -x1_8-774753035696222812246550465991661429811876498119000701158671119781263705458789269528749946767282632950719771171872558481991660348073571903818116633358984530218174409978839068897490413739637630829196939504767425837681239890499447034573005902108648219849849027980275388447363988731221559571132717943827021854747689420220378070681226929192343223127201020297746531709266956486108587145698544160593776698727548251916392870265088469574943046181966512303340557985160374122629564471402177688632705821568289084587235801620928822801/22382768364031353003183300409829275486928460080020829207538739245483860443473505625460388402468858035595463653031551373807656958676434614782746572560368117653180394077958072516360286090719405972745266136427874593301328737702802253219190855281827468054495048944859638617519754396125, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((delta*x3_7-gama_0*x2_7-7*gama_1*x2_6-21*gama_2*x2_5-35*gama_3*x2_4-35*gama_4*x2_3-21*gama_5*x2_2-7*gama_6*x2_1-gama_7*x2_0)*x4_0+(7*delta*x3_6-7*gama_0*x2_6-42*gama_1*x2_5-105*gama_2*x2_4-140*gama_3*x2_3-105*gama_4*x2_2-42*gama_5*x2_1-7*gama_6*x2_0)*x4_1+(-x2_0*x4_7-7*x2_1*x4_6-21*x2_2*x4_5-35*x2_3*x4_4-35*x2_4*x4_3-21*x2_5*x4_2)*gama_0+(x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2)*delta+(-7*gama_1*x4_6-21*gama_2*x4_5-35*gama_3*x4_4-35*gama_4*x4_3-21*gama_5*x4_2)*x2_0+(-105*gama_1*x2_4-210*gama_2*x2_3-210*gama_3*x2_2-105*gama_4*x2_1)*x4_2+(-42*x2_1*x4_5-105*x2_2*x4_4-140*x2_3*x4_3)*gama_1+(-105*gama_2*x4_4-140*gama_3*x4_3)*x2_1-210*gama_2*x2_2*x4_3)*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, gama_7, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0+x3_7, -x1_9-6768514865153563581752291225907005485744992299099128740879916009963725168078625593375311664187631820544465531445170551329496271854305509934494393575657439819919378002197860309172669583474726348440954675515233770095152294085020742898487095709794967802118242753891014868086360236279037887833113510139087719023774492821767754737630683305367806355454547478294926066497843901891397647252747545718339837432462665360872036634633231306536727959214692814848980919218589288950889646309042834077386417154480955197360548914235572676971278190381654299334571571964171043281849390301885899426062727022103842603089/382528333583179900157307187266522255443724587280053364518908031370367750580442612647317526116550839612308102381081162806200741349681146816638355382295897005867640379027505525163413302200251165995057456050453502404626670318105310076481760470138215486368960533929095223340702357503579388766809824152699303648241011735375, -gama_1, -gama_2, -gama_3, -gama_4, -gama_5, -gama_6, -gama_7, -4162804802701724193972536962888983839536616472269290331166374310404308161672079643671510278688217508273770208969541976060578722100937510180953044709360660592822510708733425313135289851216971813535868989510096449506257208531992874808404348804342202331229797501402301354367376586858622381126371052774500232384688033250570137811502204291078903976200558564214541545638402536133998249982156678772/9753253994734495628048364376039854065324129063162228730472730401004934379793607876605785197123316057779206132847596547465619494353388138348825476201391045287156726611611270867043520551516038151176580539925-x2_7, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, gama_7, x4_7, x3_7, x2_7, x1_7, gama_6, x4_6, x3_6, x2_6, x1_6, gama_5, x4_5, x3_5, x2_5, x1_5, gama_4, x4_4, x3_4, x2_4, x1_4, gama_3, x4_3, x3_3, x2_3, x1_3, gama_2, x4_2, x3_2, x2_2, x1_2, gama_1, x4_1, x3_1, x2_1, x1_1, gama_0, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [gama, x2_7] 159
# 275, [1 = 7, -1 = 2, -6*sgm*x2_2*x4_2 = 1, -3*sgm*x2_2*x4_1 = 1, delta*sgm*x3_0*x4_3 = 1, x3_2*x4_3 = 1, -7*sgm*x2_1*x4_6 = 1, b*x1_6*x4_0 = 1, x3_7 = 1, x3_3*x4_1 = 1, -7*sgm*x2_6*x4_1 = 1, x1_5*x4_1 = 1, b*c*x1_1 = 1, delta*sgm*x3_2*x4_3 = 1, x1_5*c = 1, x3_0*x4_5 = 1, b*x1_2*x4_1 = 1, x1_7*x4_2 = 1, z_aux*x3_0*c = 1, beta*x2_3 = 1, -sgm*x2_0*x4_7 = 1, x3_2 = 1, x3_6*x4_1 = 1, -x2_2 = 1, x2_3 = 1, x3_3 = 1, b*x1_2*x4_2 = 1, x3_2*x4_5 = 1, delta*x3_4 = 1, x3_2*x4_4 = 1, delta*sgm*x3_2*x4_4 = 1, x2_6 = 1, b*x1_2*x4_5 = 1, x3_0*x4_1 = 1, b*x1_5*x4_3 = 1, delta*sgm*x3_6*x4_0 = 1, b*x1_4*x4_0 = 1, x1_5*x4_3 = 1, x1_1*x4_2 = 1, x1_1*x4_8 = 1, delta*sgm*x3_4*x4_0 = 1, x1_6*x4_0 = 1, delta*sgm*x3_4*x4_3 = 1, -21*sgm*x2_5*x4_2 = 1, z_aux*x3_0*x4_0 = 1, delta*sgm*x3_7*x4_0 = 1, -x1_7 = 1, delta*x3_3 = 1, -4*sgm*x2_3*x4_1 = 1, x1_8*x4_1 = 1, b*x1_5*x4_1 = 1, -x1_6 = 1, b*x1_0*x4_2 = 1, x3_1*x4_4 = 1, -sgm*x2_6*x4_0 = 1, x3_3*x4_3 = 1, b*x1_5*x4_0 = 1, b*x1_1*x4_2 = 1, x2_2 = 1, x1_6*x4_2 = 1, -2*sgm*x2_1*x4_1 = 1, x1_1*x4_1 = 1, delta*sgm*x3_1*x4_0 = 1, -10*sgm*x2_2*x4_3 = 1, x1_4*c = 1, x2_4 = 1, -sgm*x2_0*x4_4 = 1, x3_1 = 1, delta*sgm*x3_1*x4_6 = 1, x3_1*x4_6 = 1, x3_1*x4_5 = 1, x1_4*x4_3 = 1, x3_0*x4_4 = 1, -x1_2 = 1, -35*sgm*x2_4*x4_3 = 1, x3_3*x4_4 = 1, delta*sgm*x3_0*x4_1 = 1, x1_3*x4_5 = 1, delta*x3_0 = 1, x3_6*x4_2 = 1, delta*sgm*x3_1*x4_4 = 1, -4*sgm*x2_1*x4_3 = 1, x1_1*c = 1, x1_2*x4_1 = 1, b*x1_1*x4_6 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_3*x4_1 = 1, b*x1_0*x4_5 = 1, delta*sgm*x3_2*x4_0 = 1, x1_6*c = 1, x1_5*x4_4 = 1, beta*x2_2 = 1, x3_3*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, x1_4*x4_0 = 1, x3_4*x4_2 = 1, delta*sgm*x3_5*x4_1 = 1, x3_4*x4_3 = 1, -x2_1 = 1, delta*sgm*x3_0*x4_4 = 1, x1_2*x4_0 = 1, b*x1_7*x4_1 = 1, b*x1_4*x4_3 = 1, delta*x3_5 = 1, x1_9*c = 1, -x2_4 = 1, b*c*x1_7 = 1, delta*sgm*x3_3*x4_1 = 1, -sgm*x2_1*x4_0 = 1, -x4_0*sgm = 1, -alpha*x1_3 = 1, b*x1_0*x4_1 = 1, x1_8*x4_0 = 1, x1_4*x4_2 = 1, delta*sgm*x3_3*x4_4 = 1, x1_5*x4_0 = 1, x1_1*x4_4 = 1, -3*sgm*x2_1*x4_2 = 1, x1_1*x4_3 = 1, x3_1*x4_1 = 1, x1_1*x4_0 = 1, -sgm*x2_0*x4_0 = 1, beta*x2_6 = 1, -alpha*x1_4 = 1, x3_5*x4_2 = 1, -5*sgm*x2_4*x4_1 = 1, beta*x2_0 = 1, x2_5 = 1, x3_5*x4_1 = 1, x1_5*x4_2 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_7 = 1, delta*x3_2 = 1, x2_1 = 1, -sgm*x2_0*x4_6 = 1, b*c*x1_6 = 1, -x1_4 = 1, b*c*x1_3 = 1, x3_0*x4_8 = 1, x3_7*x4_1 = 1, b*x1_4*x4_2 = 1, -5242672485478177444659387929538951251864037998431138078655813910570763801400630223669703878942695805601358660228072931893481643120033360259685914626310457154998773153/3446941865585422159032084916215751002267832245135820694126174656986797068950014823185448279 = 1, x1_6*x4_3 = 1, b*x1_7*x4_0 = 1, -sgm*x2_3*x4_0 = 1, x3_0*x4_2 = 1, b*x1_2*x4_3 = 1, b*x1_5*x4_2 = 1, delta*sgm*x3_0*x4_5 = 1, x1_2*x4_3 = 1, x1_6*x4_1 = 1, x1_3*x4_4 = 1, delta*sgm*x3_2*x4_1 = 1, b*x1_1*x4_7 = 1, delta*sgm*x3_4*x4_1 = 1, b*c*x1_0 = 1, delta*sgm*x3_0*x4_2 = 1, x1_2*c = 1, -sgm*x2_4*x4_0 = 1, x1_2*x4_4 = 1, b*x1_1*x4_1 = 1, -5*sgm*x2_1*x4_4 = 1, -15115441682291351469493999796889194256444309189702658997423903023753904476681805540674372864157150361997078285344947881840603112801189118660874555200030360184490623136736094978269849325660690618500256904101614451942602254596539800869941050913585838514524194595778754542279355296272339717990213237690915121173658700936742021099515999046296336730983160044154731834/55966407153616880985053572623591423628062510465785284303202635187864245602512659742542807287088833506264938261188496067030530498949786507070432964268747772403978351984641552701150465895475535612767 = 1, x1_2*x4_2 = 1, -21*sgm*x2_2*x4_5 = 1, x3_5 = 1, -alpha*x1_1 = 1, -285178703914060845731182048885212093993751879847955314363/8614621920582149333 = 1, -10*sgm*x2_3*x4_2 = 1, beta*x2_1 = 1, b*x1_0*x4_7 = 1, x1_3*x4_0 = 1, x3_4 = 1, -x2_6 = 1, delta*sgm*x3_0*x4_0 = 1, -x2_5 = 1, x1_2*x4_6 = 1, delta*sgm*x3_6*x4_1 = 1, x1_8*c = 1, x1_3*x4_3 = 1, delta*x3_6 = 1, x3_4*x4_1 = 1, b*c*x1_4 = 1, x1_9*x4_0 = 1, x1_3*x4_2 = 1, -15*sgm*x2_4*x4_2 = 1, x3_0*x4_6 = 1, beta*x2_5 = 1, -alpha*x1_2 = 1, x3_5*x4_3 = 1, delta*sgm*x3_4*x4_2 = 1, b*x1_0*x4_6 = 1, -35*sgm*x2_3*x4_4 = 1, delta*sgm*x3_3*x4_0 = 1, x1_4*x4_5 = 1, b*x1_4*x4_1 = 1, -x1_8 = 1, delta*x3_1 = 1, b*x1_3*x4_0 = 1, -sgm*x2_0*x4_2 = 1, -sgm*x2_5*x4_0 = 1, delta*sgm*x3_2*x4_2 = 1, -x1_9 = 1, b*c*x1_8 = 1, -x2_3 = 1, b*x1_0*x4_8 = 1, b*x1_6*x4_1 = 1, -5629317267317686165489021192098964322447737016015043776670822618707790786426718509586239677710877917116627651791041876428312678539839048764250576436518854083802900965418286907046686072707809546263116255395066494176408946276651216573074061161588942758327385018249058296052529508589053836325/153246013764087925164908048586910532790333032454155717482009559419018673359698329123527617565138776176515502005311383925632499363141958349866424164716281220390053 = 1, x1_3*x4_1 = 1, delta*sgm*x3_3*x4_2 = 1, x3_6 = 1, x1_1*x4_7 = 1, x1_7*x4_0 = 1, x3_2*x4_2 = 1, x3_1*x4_2 = 1, -x1_3 = 1, x3_2*x4_6 = 1, -x1_5 = 1, delta*sgm*x3_1*x4_1 = 1, x3_3*x4_5 = 1, -20*sgm*x2_3*x4_3 = 1, -6*sgm*x2_1*x4_5 = 1, b*x1_2*x4_4 = 1, b*x1_4*x4_4 = 1, delta*sgm*x3_1*x4_2 = 1, x1_3*x4_6 = 1, b*x1_0*x4_0 = 1, b*x1_6*x4_2 = 1, x1_1*x4_6 = 1, x3_0*x4_3 = 1, x1_4*x4_1 = 1, -sgm*x2_0*x4_1 = 1, -x1_0 = 1, -sgm*x2_0*x4_5 = 1, -alpha*x1_0 = 1, delta*sgm*x3_5*x4_2 = 1, b*x1_1*x4_3 = 1, b*x1_1*x4_0 = 1, x3_0*x4_7 = 1, -sgm*x2_0*x4_3 = 1, x3_4*x4_4 = 1, x3_1*x4_7 = 1, b*x1_8*x4_0 = 1, b*x1_0*x4_4 = 1, x1_2*x4_5 = 1, -alpha*x1_6 = 1, -6*sgm*x2_5*x4_1 = 1, -15*sgm*x2_2*x4_4 = 1, x1_2*x4_7 = 1, x1_7*c = 1, b*x1_1*x4_4 = 1, x3_1*x4_3 = 1, x3_2*x4_1 = 1, b*x1_3*x4_4 = 1, x1_3*c = 1, delta*sgm*x3_2*x4_5 = 1, b*c*x1_5 = 1, b*x1_3*x4_5 = 1, delta*sgm*x3_3*x4_3 = 1, b*x1_1*x4_5 = 1, -sgm*x2_2*x4_0 = 1, b*x1_0*x4_3 = 1, -alpha*x1_5 = 1, beta*x2_4 = 1, b*x1_2*x4_6 = 1, x1_7*x4_1 = 1, x1_1*x4_5 = 1, b*x1_2*x4_0 = 1, b*x1_3*x4_2 = 1, -x2_0 = 1, -x1_1 = 1, delta*sgm*x3_1*x4_3 = 1, b*x1_3*x4_3 = 1, delta*sgm*x3_0*x4_6 = 1, b*c*x1_2 = 1]
# 282, -5.588688399
# 2
# [gama = 43, x2_7 = 1]
# [gama = [4, 4, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2], x2_7 = [1]]
# [gama = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], x2_7 = [1, 1]]