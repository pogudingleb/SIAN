infolevel[Groebner]:=10;
Et_hat := [560729108360377062-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-22538918461492368546395288238751182382825568290974210083/6529049193443110367, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+909645178487798493143674308718378109017000286013584491387211644334270755434073396443027445652874041247426370788/42801487285296571265147318659523091383573582866546718477, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-256985731818436776329711160753050235661157401508292498561355661169944995528365096460152096331539371334418334258660739065827578219541722460839642961861220955297205386396/1964110059043862199732163782272535089673368563937053012701391385388368983797692305061637990209, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+72601369124167647533707473255778983443678533628407288073858139993612072344087956752560882003599481705195347563843782722246559575168494223765607025128687553161797078032261027622207346404900191454359944431824088492791170711202/90130707335548928192026900349087801994258968577311055426338455990743146045835473116642243440269918354120658589615620364343649757253, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5-3400332316712674199059850033232192093235419272280925965897104209848205044802100448321230558188729858973707839033023713156659328752926704854681169428899019468699963332041672467914068230824069549300630326563425438463732980833120258566843515167778603630079698688245787248443653738979035909173244/4135992465086682623154347014943156193093103837512160325038556785800537331859623044344945056979813343594152665777845698770339544069295837621839974198567789717726116986201, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6+1879697083603352133593783423169558687175292787470936940827775798696665808350408265911493118949383859342964935383459888376929453712806211096991287447544466142036652301395879280881226714987492542538464232039167094159298643934532058895215262759739068648362566112696030819910827847491354858106198800880187549965216643841788747516361947814099372959492512703471362095312/3873384571923615183961426505402747757796360051323339940641358706203957078453137254934790507331837640334065331755878181403315113250079718370052939639716990093811096055917091458623284391192147388370536441933, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7+567008147084517956783420450292510191091458048066418247251575963007071269753455211166836031173562876217266219145614772685892540896617024189746967884998887101215202194720872381909008056389548913951100332842880714398099633135916192371319700935357859130817771801718242847918527940400744325164164225687181566074113935112704192716640302158583033216419693945130419115431401656917333171433041922221163378401896053909456808153591313728593440532372/1244215530335497745482870304152891178932086544519042769077264042091644871319594990251927313235595939325927912588082631999186877128301981963228782675126732459745662464066589762720649251696124178342129951772546901811052863748269794763574274248327, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8-10471215976645838564959548673460043562491527152985911865582174194478579652154010888416064069210642064603932774355566961679979393437347538457965307414351041944964509230044706106221669505612966436105512751789763594980233772990214304840338871984416093266744462604228325496410503181860196646734034970900745215802624126208047984397730400549331487516321899410959163544920696961156279221004450520268246719449184874033023165635934202634498313768680606425637117927924451381071396480098340884917852709664876203099801954864/57095591619547430054572675768813928987212526350215594234590109116871880440247779961346889412030412067658139022149789106878985900102466958581570135491655362030141826808426333308137304260373701927015252242507767535018709400937729425426112980904931633877566596276138818619971450602659, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9-1097371831993845233385774822555959246505250784619911281064787485460941916056606951227228245875034027086175214603748516690430276377232027378284364200982257282224620670819541231942334315035054505166832601962542698208114081006641431697907836903678585502842557003032854269723314460759355314642811898555997595524608933883149091433735805705165425108698035439540680252205292464344380478839478697560904940095534229775667885941242344409635065803877539534476599636546111807886173574258282757216569054226329176581460837763764825253381066428293488727761718273641140649170625044733378932224351715392/2620049744522248443959553989699751749382930826969749977659741504803782956118352339072534514049914660146719068797601522281082896153142088066622572604435228540129238388711984822254210090880151533396881100125083271470101989383794104310643314839036778764872548348990912817310206971437401398150276721324402074900496948953903, 1137739345349326738012978135044394514136528963286422612115278148363905956139643309140993905431510785183387800201730887744748403884922299762225484967476693752703099436955196904503052433254492934409103154066948639234400816640494539914780829501462/90130707335548928192026900349087801994258968577311055426338455990743146045835473116642243440269918354120658589615620364343649757253-x2_5, -116263644511453343538494329649474381451940920314397227545814312289748815418719139933424873421194560374842385916401088088038559465216175175261530726339334360648129695248398800034684455493616094914759610502430480292109055366557369876973770504781366536913486674443295481939689669844682163258237164885557134582997867109018862438179432/4135992465086682623154347014943156193093103837512160325038556785800537331859623044344945056979813343594152665777845698770339544069295837621839974198567789717726116986201-x3_7, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [x2_5, x3_7] 24
# 274, [1 = 8, -1 = 2, delta*sgm*x3_0*x4_3 = 1, x3_2*x4_3 = 1, x4_1 = 1, b*x1_6*x4_0 = 1, x3_3*x4_1 = 1, x1_5*x4_1 = 1, b*c*x1_1 = 1, -15*gama*sgm*x2_2*x4_4 = 1, delta*sgm*x3_2*x4_3 = 1, x1_5*c = 1, x3_0*x4_5 = 1, b*x1_2*x4_1 = 1, x1_7*x4_2 = 1, z_aux*x3_0*c = 1, beta*x2_3 = 1, x3_2 = 1, x3_6*x4_1 = 1, x2_3 = 1, x3_3 = 1, b*x1_2*x4_2 = 1, -gama*x2_4 = 1, x3_2*x4_5 = 1, delta*x3_4 = 1, x3_2*x4_4 = 1, delta*sgm*x3_2*x4_4 = 1, x2_6 = 1, b*x1_2*x4_5 = 1, -gama*sgm*x2_0*x4_1 = 1, x3_0*x4_1 = 1, b*x1_5*x4_3 = 1, delta*sgm*x3_6*x4_0 = 1, b*x1_4*x4_0 = 1, x1_5*x4_3 = 1, x1_1*x4_2 = 1, x1_1*x4_8 = 1, delta*sgm*x3_4*x4_0 = 1, x1_6*x4_0 = 1, x2_7 = 1, -gama*sgm*x2_3*x4_0 = 1, -x4_0*gama*sgm = 1, delta*sgm*x3_4*x4_3 = 1, z_aux*x3_0*x4_0 = 1, -x1_7 = 1, delta*x3_3 = 1, x1_8*x4_1 = 1, b*x1_5*x4_1 = 1, -x1_6 = 1, -6*gama*sgm*x2_1*x4_5 = 1, b*x1_0*x4_2 = 1, x3_1*x4_4 = 1, -gama*sgm*x2_0*x4_6 = 1, x3_3*x4_3 = 1, b*x1_5*x4_0 = 1, b*x1_1*x4_2 = 1, -3*gama*sgm*x2_2*x4_1 = 1, x2_2 = 1, x1_6*x4_2 = 1, x1_1*x4_1 = 1, delta*sgm*x3_1*x4_0 = 1, -gama*sgm*x2_0*x4_4 = 1, x1_4*c = 1, x2_4 = 1, -gama*x2_3 = 1, x3_1 = 1, delta*sgm*x3_1*x4_6 = 1, x3_1*x4_6 = 1, x3_1*x4_5 = 1, x1_4*x4_3 = 1, x3_0*x4_4 = 1, -x1_2 = 1, x3_3*x4_4 = 1, delta*sgm*x3_0*x4_1 = 1, x1_3*x4_5 = 1, -3*gama*sgm*x2_1*x4_2 = 1, delta*x3_0 = 1, x3_6*x4_2 = 1, delta*sgm*x3_1*x4_4 = 1, x1_1*c = 1, x1_2*x4_1 = 1, b*x1_1*x4_6 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_3*x4_1 = 1, b*x1_0*x4_5 = 1, delta*sgm*x3_2*x4_0 = 1, x1_6*c = 1, x1_5*x4_4 = 1, beta*x2_2 = 1, x3_3*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, x1_4*x4_0 = 1, -gama*sgm*x2_6*x4_0 = 1, x3_4*x4_2 = 1, delta*sgm*x3_5*x4_1 = 1, x3_4*x4_3 = 1, delta*sgm*x3_0*x4_4 = 1, x1_2*x4_0 = 1, b*x1_7*x4_1 = 1, b*x1_4*x4_3 = 1, delta*x3_5 = 1, -21*gama*x4_2*sgm = 1, x1_9*c = 1, -gama*sgm*x2_2*x4_0 = 1, -15*gama*sgm*x2_4*x4_2 = 1, b*c*x1_7 = 1, delta*sgm*x3_3*x4_1 = 1, -alpha*x1_3 = 1, b*x1_0*x4_1 = 1, x1_8*x4_0 = 1, x1_4*x4_2 = 1, delta*sgm*x3_3*x4_4 = 1, x1_5*x4_0 = 1, x1_1*x4_4 = 1, x1_1*x4_3 = 1, x3_1*x4_1 = 1, x1_1*x4_0 = 1, beta*x2_6 = 1, -alpha*x1_4 = 1, x3_5*x4_2 = 1, -gama*sgm*x2_0*x4_3 = 1, beta*x2_0 = 1, -gama*sgm*x2_0*x4_2 = 1, x3_5*x4_1 = 1, x1_5*x4_2 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_7 = 1, delta*x3_2 = 1, x2_1 = 1, -gama*sgm*x2_0*x4_0 = 1, -20*gama*sgm*x2_3*x4_3 = 1, -2*gama*sgm*x2_1*x4_1 = 1, b*c*x1_6 = 1, -x1_4 = 1, b*c*x1_3 = 1, x3_0*x4_8 = 1, b*x1_4*x4_2 = 1, -5242672485478177444659387929538951251864037998431138078655813910570763801400630223669703878942695805601358660228072931893481643120033360259685914626310457154998773153/3446941865585422159032084916215751002267832245135820694126174656986797068950014823185448279 = 1, x1_6*x4_3 = 1, b*x1_7*x4_0 = 1, x3_0*x4_2 = 1, -7*gama*sgm*x2_6*x4_1 = 1, b*x1_2*x4_3 = 1, b*x1_5*x4_2 = 1, delta*sgm*x3_0*x4_5 = 1, x1_2*x4_3 = 1, x1_6*x4_1 = 1, -gama*sgm*x2_4*x4_0 = 1, x1_3*x4_4 = 1, delta*sgm*x3_2*x4_1 = 1, b*x1_1*x4_7 = 1, delta*sgm*x3_4*x4_1 = 1, b*c*x1_0 = 1, delta*sgm*x3_0*x4_2 = 1, -gama*x2_1 = 1, x1_2*c = 1, -gama*sgm*x2_0*x4_7 = 1, x1_2*x4_4 = 1, b*x1_1*x4_1 = 1, -15115441682291351469493999796889194256444309189702658997423903023753904476681805540674372864157150361997078285344947881840603112801189118660874555200030360184490623136736094978269849325660690618500256904101614451942602254596539800869941050913585838514524194595778754542279355296272339717990213237690915121173658700936742021099515999046296336730983160044154731834/55966407153616880985053572623591423628062510465785284303202635187864245602512659742542807287088833506264938261188496067030530498949786507070432964268747772403978351984641552701150465895475535612767 = 1, x1_2*x4_2 = 1, x3_5 = 1, -alpha*x1_1 = 1, -285178703914060845731182048885212093993751879847955314363/8614621920582149333 = 1, -4*gama*sgm*x2_3*x4_1 = 1, beta*x2_1 = 1, b*x1_0*x4_7 = 1, x1_3*x4_0 = 1, x3_4 = 1, delta*sgm*x3_0*x4_0 = 1, x1_2*x4_6 = 1, delta*sgm*x3_6*x4_1 = 1, x1_8*c = 1, x1_3*x4_3 = 1, delta*x3_6 = 1, x3_4*x4_1 = 1, b*c*x1_4 = 1, x1_9*x4_0 = 1, -gama = 1, beta = 1, x1_3*x4_2 = 1, -gama*sgm*x2_7*x4_0 = 1, x3_0*x4_6 = 1, -alpha*x1_2 = 1, x3_5*x4_3 = 1, delta*sgm*x3_4*x4_2 = 1, b*x1_0*x4_6 = 1, delta*sgm*x3_3*x4_0 = 1, x1_4*x4_5 = 1, b*x1_4*x4_1 = 1, -x1_8 = 1, delta*x3_1 = 1, b*x1_3*x4_0 = 1, delta*sgm*x3_2*x4_2 = 1, -gama*x2_0 = 1, -x1_9 = 1, b*c*x1_8 = 1, b*x1_0*x4_8 = 1, b*x1_6*x4_1 = 1, -5629317267317686165489021192098964322447737016015043776670822618707790786426718509586239677710877917116627651791041876428312678539839048764250576436518854083802900965418286907046686072707809546263116255395066494176408946276651216573074061161588942758327385018249058296052529508589053836325/153246013764087925164908048586910532790333032454155717482009559419018673359698329123527617565138776176515502005311383925632499363141958349866424164716281220390053 = 1, -35*gama*sgm*x2_3*x4_4 = 1, x1_3*x4_1 = 1, -4*gama*sgm*x2_1*x4_3 = 1, delta*sgm*x3_3*x4_2 = 1, -gama*x2_6 = 1, x3_6 = 1, x1_1*x4_7 = 1, x1_7*x4_0 = 1, x3_2*x4_2 = 1, x3_1*x4_2 = 1, -x1_3 = 1, -10*gama*sgm*x2_3*x4_2 = 1, x3_2*x4_6 = 1, -x1_5 = 1, -gama*x2_2 = 1, delta*sgm*x3_1*x4_1 = 1, x3_3*x4_5 = 1, -10*gama*sgm*x2_2*x4_3 = 1, b*x1_2*x4_4 = 1, b*x1_4*x4_4 = 1, delta*sgm*x3_1*x4_2 = 1, -21*gama*sgm*x2_2*x4_5 = 1, x1_3*x4_6 = 1, b*x1_0*x4_0 = 1, b*x1_6*x4_2 = 1, x1_1*x4_6 = 1, -35*gama*sgm*x2_4*x4_3 = 1, x3_0*x4_3 = 1, x1_4*x4_1 = 1, -x1_0 = 1, -alpha*x1_0 = 1, delta*sgm*x3_5*x4_2 = 1, -gama*sgm*x2_1*x4_0 = 1, b*x1_1*x4_3 = 1, b*x1_1*x4_0 = 1, x3_0*x4_7 = 1, -5*gama*sgm*x2_1*x4_4 = 1, x3_4*x4_4 = 1, delta*x4_0*sgm = 1, x3_1*x4_7 = 1, b*x1_8*x4_0 = 1, -5*gama*sgm*x2_4*x4_1 = 1, b*x1_0*x4_4 = 1, x1_2*x4_5 = 1, -alpha*x1_6 = 1, -6*gama*x4_1*sgm = 1, x1_2*x4_7 = 1, x1_7*c = 1, b*x1_1*x4_4 = 1, x3_1*x4_3 = 1, x3_2*x4_1 = 1, -7*gama*sgm*x2_1*x4_6 = 1, b*x1_3*x4_4 = 1, x1_3*c = 1, delta*sgm*x3_2*x4_5 = 1, b*c*x1_5 = 1, b*x1_3*x4_5 = 1, delta*sgm*x3_3*x4_3 = 1, b*x1_1*x4_5 = 1, b*x1_0*x4_3 = 1, -alpha*x1_5 = 1, beta*x2_4 = 1, b*x1_2*x4_6 = 1, x1_7*x4_1 = 1, -gama*sgm*x2_0*x4_5 = 1, x1_1*x4_5 = 1, b*x1_2*x4_0 = 1, b*x1_3*x4_2 = 1, -x1_1 = 1, delta*sgm*x3_1*x4_3 = 1, b*x1_3*x4_3 = 1, -6*gama*sgm*x2_2*x4_2 = 1, delta*sgm*x3_0*x4_6 = 1, b*c*x1_2 = 1]
# 282, -5.577999884
# 2
# [x2_5 = 6, x3_7 = 3]
# [x2_5 = [4, 1, 4, 2, 2, 4], x3_7 = [4, 2, 1]]
# [x2_5 = [1, 1, 1, 1, 1, 1], x3_7 = [1, 1, 1]]