infolevel[Groebner]:=10;
Et_hat := [5737378270086194650-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-414218283699368957906113474885915622037799064716762413299/14428495344467101826, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+323529441178552851461351795422950460096102781167729799992570428034045050281877277269613281935212434636399356829/2252220200692740236640690589725450031761038718224018506, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-252695990081584491619449522117867564776894433412116505889320043620198923535712974402405299964084021177585985446506720467403200128852472872639537151962713395835257697/351560970933368888200575939397140172980943350022997003582455220904378485198423086315777586, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+50186959316847442818406983596611425607722092935295419573654019548219691931455944586985482886919998434633521446708329906252522217193854720214423711242489300874858598482765453811064666729970118083401954079047623792099904525/54877012578786714371199886418022294433453213887846965609607088017143752414519601110816422510808832424005507977599027105211066, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5-2796352915189277331301691036638247923055184980324048299354040364255596214345473245995597102851227045329724122125630539137441910907347461228682806400807680123892217389785252165995423337221378527439042569684067273971630636591782695929719510134062847480429230806877245032273065606208177097613814289/8566043328350803382159081996834816521991604353129434478698621855196154677004417354189852426670875407428519460456511901151526004079982476239421573787386280910946, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6-252727670723567056965563942214586233263415954420214932148034442028356767571257573167001460592595310187701082185084947070542605659164971834391455618278308031505704066649736502295789958209869693989114804055131998798741505560188662493225255422252282464527131456500727967970077185144304381918849119112205911200226193708775548777819065156376689301789769234266850553610648371/1337119038647304164895853337063582678918911610187322088499946743336917171620337726096484199085459768987642522736233952495608866026170239782509981550135261273180943283484966624422939900366374029226, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7+35543019728809685402478307933071597644037684993409571285270483315511623231553115235208475482122078468475124980510371460545035906790340166463672745277864710146338110873328045046212743716465224236010031268287005025591556928680995824413273220486413803207643046835731781499279841250751504973432430375640774401375396661062452585013410093115536596436717249048282506687383885473465331220196831181651006746153262319157624738485718876407267346168269423/208717987404496106821226475576407693741365498004991109128354138874963869267866164574329335907754945874985376763221236548688654627733247551746118380695032049002143421211682312776462183510148561737657838785524258782874682927699077906, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8+2209425309517045773792033196285607959625231152453955369944404737512178926665149634700291482503176164778071919862988059098115843974286494702575573498131950526446838726295933234248581788123880418953602696412095675953217289938945597238300880567975715617494637214356241159919636746867322570078741401444516114886206782788349162540325526854872004656406581881243105762841319319687069263090047387425339168952954868767732319515884474591696563723065631893903048176051403481534881685388369553248453880452523496945671449495143405/32579895287598390759140918258083966651310746687465192561487947086799672277984734490548065086116468476415450455231894021989391132619278129927902299232495422205726326360816178262525214224068938936709368671153549118264884846042218919352936253696758446744405469038728986, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9-749462506202842058759995622062039857471312609952386419692997632250030888112706779126468626291465906057923356214742070020989780518057145674498083817822240735181566114465280191791334699335803616160261894269891125585662439257791077997971523219234200558934894220848020093096671939433067249273856602060983654242509960542128008399097649970748609652424825850420738674467871906116786085603735542388243913204560430798968502129200147361997752292708474061916741655627543242166914241288453647827264267318435272988179674984635699644226803141679510564441509822678997781547663331417747544950502514953246129/5085568283550872197931662720476346034554696918643380004063496376098787175052194730430460333499673647687316232962216264014889282329349538119590862589730496566100855614241739044044054702893471814204839775500969217433617333632578091210265332700698652505142644932612139921182739017398206141053213992614466, -1173866329118750018696318964188420120541297722362802237343809516634148944487/7214247672233550913-x2_2, -42870209692830864920452262690139615963092018709460782302130192238920969668990696834740279575248249838037748803651053032110482855561171872615654107700484271480667538395305734938500163035256538576485769963267971526641319286900612860664135027647194951068372452598455258379297893626788964790083690960078026853000957264483944347639827222/4283021664175401691079540998417408260995802176564717239349310927598077338502208677094926213335437703714259730228255950575763002039991238119710786893693140455473-x3_7, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [x2_2, x3_7] 36
# 274, [1 = 8, -1 = 2, delta*sgm*x3_0*x4_3 = 1, x3_2*x4_3 = 1, x4_1 = 1, b*x1_6*x4_0 = 1, x3_3*x4_1 = 1, x1_5*x4_1 = 1, b*c*x1_1 = 1, delta*sgm*x3_2*x4_3 = 1, x1_5*c = 1, x3_0*x4_5 = 1, b*x1_2*x4_1 = 1, x1_7*x4_2 = 1, z_aux*x3_0*c = 1, beta*x2_3 = 1, x3_2 = 1, x3_6*x4_1 = 1, x2_3 = 1, x3_3 = 1, b*x1_2*x4_2 = 1, -gama*x2_4 = 1, x3_2*x4_5 = 1, delta*x3_4 = 1, x3_2*x4_4 = 1, delta*sgm*x3_2*x4_4 = 1, x2_6 = 1, b*x1_2*x4_5 = 1, -gama*sgm*x2_0*x4_1 = 1, x3_0*x4_1 = 1, b*x1_5*x4_3 = 1, delta*sgm*x3_6*x4_0 = 1, -6*gama*x4_2*sgm = 1, b*x1_4*x4_0 = 1, x1_5*x4_3 = 1, x1_1*x4_2 = 1, x1_1*x4_8 = 1, delta*sgm*x3_4*x4_0 = 1, x1_6*x4_0 = 1, x2_7 = 1, -gama*sgm*x2_3*x4_0 = 1, -x4_0*gama*sgm = 1, delta*sgm*x3_4*x4_3 = 1, z_aux*x3_0*x4_0 = 1, -x1_7 = 1, delta*x3_3 = 1, x1_8*x4_1 = 1, b*x1_5*x4_1 = 1, -x1_6 = 1, -6*gama*sgm*x2_1*x4_5 = 1, b*x1_0*x4_2 = 1, x3_1*x4_4 = 1, -gama*sgm*x2_0*x4_6 = 1, x3_3*x4_3 = 1, b*x1_5*x4_0 = 1, b*x1_1*x4_2 = 1, x1_6*x4_2 = 1, x1_1*x4_1 = 1, delta*sgm*x3_1*x4_0 = 1, -gama*sgm*x2_0*x4_4 = 1, x1_4*c = 1, x2_4 = 1, -gama*x2_3 = 1, x3_1 = 1, delta*sgm*x3_1*x4_6 = 1, x3_1*x4_6 = 1, x3_1*x4_5 = 1, x1_4*x4_3 = 1, x3_0*x4_4 = 1, -x1_2 = 1, x3_3*x4_4 = 1, delta*sgm*x3_0*x4_1 = 1, x1_3*x4_5 = 1, -3*gama*sgm*x2_1*x4_2 = 1, delta*x3_0 = 1, x3_6*x4_2 = 1, -6*gama*sgm*x2_5*x4_1 = 1, delta*sgm*x3_1*x4_4 = 1, x1_1*c = 1, x1_2*x4_1 = 1, b*x1_1*x4_6 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_3*x4_1 = 1, b*x1_0*x4_5 = 1, delta*sgm*x3_2*x4_0 = 1, x1_6*c = 1, x1_5*x4_4 = 1, x3_3*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, x1_4*x4_0 = 1, -gama*sgm*x2_6*x4_0 = 1, x3_4*x4_2 = 1, delta*sgm*x3_5*x4_1 = 1, x3_4*x4_3 = 1, delta*sgm*x3_0*x4_4 = 1, x1_2*x4_0 = 1, b*x1_7*x4_1 = 1, b*x1_4*x4_3 = 1, delta*x3_5 = 1, x1_9*c = 1, -15*gama*sgm*x2_4*x4_2 = 1, b*c*x1_7 = 1, delta*sgm*x3_3*x4_1 = 1, -alpha*x1_3 = 1, b*x1_0*x4_1 = 1, x1_8*x4_0 = 1, x1_4*x4_2 = 1, delta*sgm*x3_3*x4_4 = 1, x1_5*x4_0 = 1, x1_1*x4_4 = 1, x1_1*x4_3 = 1, x3_1*x4_1 = 1, x1_1*x4_0 = 1, beta*x2_6 = 1, -alpha*x1_4 = 1, x3_5*x4_2 = 1, -gama*sgm*x2_0*x4_3 = 1, beta*x2_0 = 1, x2_5 = 1, -gama*sgm*x2_0*x4_2 = 1, x3_5*x4_1 = 1, x1_5*x4_2 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_7 = 1, delta*x3_2 = 1, x2_1 = 1, -gama*sgm*x2_0*x4_0 = 1, -20*gama*sgm*x2_3*x4_3 = 1, -2*gama*sgm*x2_1*x4_1 = 1, b*c*x1_6 = 1, -15*gama*x4_4*sgm = 1, -x1_4 = 1, b*c*x1_3 = 1, x3_0*x4_8 = 1, b*x1_4*x4_2 = 1, -5242672485478177444659387929538951251864037998431138078655813910570763801400630223669703878942695805601358660228072931893481643120033360259685914626310457154998773153/3446941865585422159032084916215751002267832245135820694126174656986797068950014823185448279 = 1, x1_6*x4_3 = 1, b*x1_7*x4_0 = 1, x3_0*x4_2 = 1, -7*gama*sgm*x2_6*x4_1 = 1, b*x1_2*x4_3 = 1, b*x1_5*x4_2 = 1, delta*sgm*x3_0*x4_5 = 1, x1_2*x4_3 = 1, x1_6*x4_1 = 1, -gama*sgm*x2_4*x4_0 = 1, x1_3*x4_4 = 1, delta*sgm*x3_2*x4_1 = 1, b*x1_1*x4_7 = 1, delta*sgm*x3_4*x4_1 = 1, b*c*x1_0 = 1, delta*sgm*x3_0*x4_2 = 1, -gama*x2_1 = 1, x1_2*c = 1, -gama*sgm*x2_0*x4_7 = 1, x1_2*x4_4 = 1, b*x1_1*x4_1 = 1, -15115441682291351469493999796889194256444309189702658997423903023753904476681805540674372864157150361997078285344947881840603112801189118660874555200030360184490623136736094978269849325660690618500256904101614451942602254596539800869941050913585838514524194595778754542279355296272339717990213237690915121173658700936742021099515999046296336730983160044154731834/55966407153616880985053572623591423628062510465785284303202635187864245602512659742542807287088833506264938261188496067030530498949786507070432964268747772403978351984641552701150465895475535612767 = 1, x1_2*x4_2 = 1, x3_5 = 1, -alpha*x1_1 = 1, -285178703914060845731182048885212093993751879847955314363/8614621920582149333 = 1, -10*gama*x4_3*sgm = 1, -4*gama*sgm*x2_3*x4_1 = 1, beta*x2_1 = 1, b*x1_0*x4_7 = 1, x1_3*x4_0 = 1, x3_4 = 1, delta*sgm*x3_0*x4_0 = 1, x1_2*x4_6 = 1, delta*sgm*x3_6*x4_1 = 1, x1_8*c = 1, x1_3*x4_3 = 1, delta*x3_6 = 1, -21*gama*x4_5*sgm = 1, x3_4*x4_1 = 1, b*c*x1_4 = 1, x1_9*x4_0 = 1, -gama = 1, beta = 1, x1_3*x4_2 = 1, -gama*sgm*x2_7*x4_0 = 1, x3_0*x4_6 = 1, beta*x2_5 = 1, -alpha*x1_2 = 1, -3*gama*x4_1*sgm = 1, x3_5*x4_3 = 1, delta*sgm*x3_4*x4_2 = 1, b*x1_0*x4_6 = 1, delta*sgm*x3_3*x4_0 = 1, x1_4*x4_5 = 1, b*x1_4*x4_1 = 1, -x1_8 = 1, delta*x3_1 = 1, b*x1_3*x4_0 = 1, delta*sgm*x3_2*x4_2 = 1, -gama*x2_0 = 1, -x1_9 = 1, b*c*x1_8 = 1, b*x1_0*x4_8 = 1, b*x1_6*x4_1 = 1, -5629317267317686165489021192098964322447737016015043776670822618707790786426718509586239677710877917116627651791041876428312678539839048764250576436518854083802900965418286907046686072707809546263116255395066494176408946276651216573074061161588942758327385018249058296052529508589053836325/153246013764087925164908048586910532790333032454155717482009559419018673359698329123527617565138776176515502005311383925632499363141958349866424164716281220390053 = 1, -35*gama*sgm*x2_3*x4_4 = 1, x1_3*x4_1 = 1, -4*gama*sgm*x2_1*x4_3 = 1, delta*sgm*x3_3*x4_2 = 1, -gama*x2_6 = 1, x3_6 = 1, x1_1*x4_7 = 1, x1_7*x4_0 = 1, x3_2*x4_2 = 1, x3_1*x4_2 = 1, -x1_3 = 1, -10*gama*sgm*x2_3*x4_2 = 1, x3_2*x4_6 = 1, -x1_5 = 1, delta*sgm*x3_1*x4_1 = 1, x3_3*x4_5 = 1, b*x1_2*x4_4 = 1, b*x1_4*x4_4 = 1, delta*sgm*x3_1*x4_2 = 1, x1_3*x4_6 = 1, b*x1_0*x4_0 = 1, b*x1_6*x4_2 = 1, x1_1*x4_6 = 1, -35*gama*sgm*x2_4*x4_3 = 1, x3_0*x4_3 = 1, x1_4*x4_1 = 1, -x1_0 = 1, -alpha*x1_0 = 1, delta*sgm*x3_5*x4_2 = 1, -gama*sgm*x2_1*x4_0 = 1, b*x1_1*x4_3 = 1, b*x1_1*x4_0 = 1, x3_0*x4_7 = 1, -5*gama*sgm*x2_1*x4_4 = 1, x3_4*x4_4 = 1, delta*x4_0*sgm = 1, x3_1*x4_7 = 1, b*x1_8*x4_0 = 1, -5*gama*sgm*x2_4*x4_1 = 1, b*x1_0*x4_4 = 1, x1_2*x4_5 = 1, -alpha*x1_6 = 1, x1_2*x4_7 = 1, x1_7*c = 1, b*x1_1*x4_4 = 1, x3_1*x4_3 = 1, x3_2*x4_1 = 1, -7*gama*sgm*x2_1*x4_6 = 1, b*x1_3*x4_4 = 1, x1_3*c = 1, delta*sgm*x3_2*x4_5 = 1, -gama*x2_5 = 1, b*c*x1_5 = 1, b*x1_3*x4_5 = 1, delta*sgm*x3_3*x4_3 = 1, b*x1_1*x4_5 = 1, b*x1_0*x4_3 = 1, -alpha*x1_5 = 1, beta*x2_4 = 1, -21*gama*sgm*x2_5*x4_2 = 1, b*x1_2*x4_6 = 1, x1_7*x4_1 = 1, -gama*sgm*x2_0*x4_5 = 1, x1_1*x4_5 = 1, b*x1_2*x4_0 = 1, b*x1_3*x4_2 = 1, -x1_1 = 1, delta*sgm*x3_1*x4_3 = 1, b*x1_3*x4_3 = 1, delta*sgm*x3_0*x4_6 = 1, b*c*x1_2 = 1, -gama*sgm*x2_5*x4_0 = 1]
# 282, -5.577999884
# 2
# [x2_2 = 9, x3_7 = 3]
# [x2_2 = [4, 1, 4, 2, 2, 4, 4, 4, 4], x3_7 = [4, 2, 1]]
# [x2_2 = [1, 1, 1, 1, 1, 1, 1, 1, 1], x3_7 = [1, 1, 1]]