infolevel[Groebner]:=10;
Et_hat := [7509156893425993421-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-121932436839470312242007708578463646597584880472058984993/8891677261256037051, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+345439778776545236380270465251106173448306782112791832320049070545716911480287602683141544860347832957178184878/13794069197847723732206765680398222449082572809124990683, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-14679683747929623847859234818878380624248819428168541772616858694860360939274025069215194001039939704492911428836533678228457830825358116969262500556683298594785313036/320990640085602255634954747263831507832903912769370090670757897085345325410529499975735548085, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+3232797375681321167566346066358881167896522875107431982022537494374679679290759499905248721204426122591316565261284461696382503525074404252932570545404621531729121959103044074180929033934009207098015340078865851099947062104/7469513857349729915656050900413090785534515150602154369366759418398390696511599656959352746282387392135736403642884859210550502075, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5+842595605370328085228855197689145058685978484809088141157503599208783659363408410034998558000575113535228050851827119890214325857799745918766567316149950533923163664414346764366693019577150216751907184970781443282254399516277532789049218371283052850626711119149355678865652439409020283646833747824/34763404472018568622279842816892801390561988093374000158797781555111972817611354092650062445932247684650351188678629929554532115798096972890569616477918734805265505425, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6-2063100316681635981581320417688139194807333318847921221123476062541487365411447552455983551349799759704520269934096495656624548462298813048132012029780938069160772020956094843177437371332589454730422737309876526086277613970486329194735057606025388694407918759064581025030798624790006581120861753941882024367371794399083656420785710520801144703867967049306052038287230534016/808951100141574940405875841274128064245103917781793324602551728136081997794847770147920417127093864467440440409388286171701615491282384580409232698335513400459123658436203506850348020052843388401496185375, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7+1377344062024984278237050946347027788875935114085108284727310221834419598046093309273304557509597991866521545227929490140882603799547292436301835819835365960666188088877259439444745307446481584386321718672458077662135193443278323801635481254328471744133281453330237899921400969768609278954298813141833115244083571212666224922567081763571477134050780164180971915002202241851679192915082525024341934747620764888464087857392709853262105718012022293792/18824447500445458197996378831630768336529590955552597003233964260372800408392468637460804308328392176833193056185501335117689483969932993772374800952185652410699042087587802567831982765210429761456871262359036127921640948593535564220756760625, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8+4504224244903037175768750955945456802039645802306238447616500498543084977610386865558048179701366001898752515948095142381831023140870156984238517544577493131592617986220174019799504618159700414055128244475186288338660131033254143936286443463657917210455346337138733118207690696438880496048991814737675375485437323763510538373972218658043713465395722158944782503969649386925143328886325093609664454244093188709198932085873424532739342322035395740451661227665780447382774618523154080468333667867664428534863432689258883297184/438048509526732292337971293952873715305419477979832803490648377489456313340993560904646706887826581937673597700119981051566791568815472724015225116429827951668436727865667676291971333833547484455888315277528407812542466224256351163412021101303221130237696457442433827346371509375, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9-16461420935576005449211555161843272786559718291309910950838213847437753118540472727625283561569274342020166061933877967162662014817887269758167326524718410379690691927273408828722276496716919536976820466469881985022536384049249917290420452283618339492347364809088623263440034252051265976304468173154038284318163461067804110351527954875285187003527374360302487018249828263548780890465338407581730825732814439753658169636420836293334088925447228989292671177820465872310905408516633839169784194644411427068172276018645210054706426978604508065792213957011708388613188481562394672136927971455319576198176/10193472966154831677276035462249455387502914787213036744466481577411236695793324454189652853842815400454071977292261644838809155441899984498561859056257938548819648426421388070975083579622196555592974911809957659914328057587661421627011849371645339695355668665484948835894798453111063342863039000451894954982522140625, -2176678148787437124141112033142203195431120196963549181366725124065642083867316103266373915391995814065503609120490608807320217954634730360963971541264362615561190263094158467128993040336/962971920256806766904864241791494523498711738308110272012273691256035976231588499927206644255-x2_4, 13655545985959986983790793640276760558259936422369231313301953100839280098402137638411476381520639635401254345979078437990186282293682996094691681415663043526211007615812842488429473927184531575551351946157397312814642402014165232093862211040315817053714399735192/7469513857349729915656050900413090785534515150602154369366759418398390696511599656959352746282387392135736403642884859210550502075-x3_6, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [x2_4, x3_6] 36
# 274, [1 = 8, -1 = 2, delta*sgm*x3_0*x4_3 = 1, x3_2*x4_3 = 1, x4_1 = 1, b*x1_6*x4_0 = 1, x3_7 = 1, x3_3*x4_1 = 1, x1_5*x4_1 = 1, b*c*x1_1 = 1, -15*gama*sgm*x2_2*x4_4 = 1, delta*sgm*x3_2*x4_3 = 1, x1_5*c = 1, x3_0*x4_5 = 1, b*x1_2*x4_1 = 1, x1_7*x4_2 = 1, z_aux*x3_0*c = 1, beta*x2_3 = 1, x3_2 = 1, x2_3 = 1, x3_3 = 1, b*x1_2*x4_2 = 1, x3_2*x4_5 = 1, delta*x3_4 = 1, x3_2*x4_4 = 1, delta*sgm*x3_2*x4_4 = 1, -5*gama*x4_1*sgm = 1, delta*x4_1*sgm = 1, x2_6 = 1, b*x1_2*x4_5 = 1, -gama*sgm*x2_0*x4_1 = 1, x3_0*x4_1 = 1, b*x1_5*x4_3 = 1, b*x1_4*x4_0 = 1, x1_5*x4_3 = 1, x1_1*x4_2 = 1, x1_1*x4_8 = 1, delta*sgm*x3_4*x4_0 = 1, x1_6*x4_0 = 1, x2_7 = 1, -gama*sgm*x2_3*x4_0 = 1, -x4_0*gama*sgm = 1, delta*sgm*x3_4*x4_3 = 1, z_aux*x3_0*x4_0 = 1, delta*sgm*x3_7*x4_0 = 1, -x1_7 = 1, delta*x3_3 = 1, x1_8*x4_1 = 1, b*x1_5*x4_1 = 1, -x1_6 = 1, -6*gama*sgm*x2_1*x4_5 = 1, b*x1_0*x4_2 = 1, x3_1*x4_4 = 1, -gama*sgm*x2_0*x4_6 = 1, x3_3*x4_3 = 1, b*x1_5*x4_0 = 1, b*x1_1*x4_2 = 1, -3*gama*sgm*x2_2*x4_1 = 1, x2_2 = 1, x1_6*x4_2 = 1, x1_1*x4_1 = 1, delta*sgm*x3_1*x4_0 = 1, -gama*sgm*x2_0*x4_4 = 1, x1_4*c = 1, -gama*x2_3 = 1, x3_1 = 1, -35*gama*x4_3*sgm = 1, delta*sgm*x3_1*x4_6 = 1, x3_1*x4_6 = 1, x3_1*x4_5 = 1, x1_4*x4_3 = 1, x3_0*x4_4 = 1, -x1_2 = 1, x3_3*x4_4 = 1, delta*sgm*x3_0*x4_1 = 1, x1_3*x4_5 = 1, -3*gama*sgm*x2_1*x4_2 = 1, delta*x3_0 = 1, -6*gama*sgm*x2_5*x4_1 = 1, delta*sgm*x3_1*x4_4 = 1, x1_1*c = 1, x1_2*x4_1 = 1, b*x1_1*x4_6 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_3*x4_1 = 1, b*x1_0*x4_5 = 1, delta*sgm*x3_2*x4_0 = 1, x1_6*c = 1, x1_5*x4_4 = 1, beta*x2_2 = 1, x3_3*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, x1_4*x4_0 = 1, -gama*sgm*x2_6*x4_0 = 1, x4_2 = 1, x3_4*x4_2 = 1, delta*sgm*x3_5*x4_1 = 1, x3_4*x4_3 = 1, delta*sgm*x3_0*x4_4 = 1, x1_2*x4_0 = 1, b*x1_7*x4_1 = 1, b*x1_4*x4_3 = 1, delta*x3_5 = 1, x1_9*c = 1, -gama*sgm*x2_2*x4_0 = 1, b*c*x1_7 = 1, delta*sgm*x3_3*x4_1 = 1, -alpha*x1_3 = 1, b*x1_0*x4_1 = 1, x1_8*x4_0 = 1, x1_4*x4_2 = 1, delta*sgm*x3_3*x4_4 = 1, x1_5*x4_0 = 1, x1_1*x4_4 = 1, x1_1*x4_3 = 1, x3_1*x4_1 = 1, x1_1*x4_0 = 1, beta*x2_6 = 1, -alpha*x1_4 = 1, x3_5*x4_2 = 1, -gama*sgm*x2_0*x4_3 = 1, beta*x2_0 = 1, x2_5 = 1, -gama*sgm*x2_0*x4_2 = 1, x3_5*x4_1 = 1, x1_5*x4_2 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_7 = 1, delta*x3_2 = 1, x2_1 = 1, -gama*sgm*x2_0*x4_0 = 1, -20*gama*sgm*x2_3*x4_3 = 1, -2*gama*sgm*x2_1*x4_1 = 1, b*c*x1_6 = 1, -x1_4 = 1, b*c*x1_3 = 1, x3_0*x4_8 = 1, x3_7*x4_1 = 1, b*x1_4*x4_2 = 1, -5242672485478177444659387929538951251864037998431138078655813910570763801400630223669703878942695805601358660228072931893481643120033360259685914626310457154998773153/3446941865585422159032084916215751002267832245135820694126174656986797068950014823185448279 = 1, x1_6*x4_3 = 1, b*x1_7*x4_0 = 1, x3_0*x4_2 = 1, -7*gama*sgm*x2_6*x4_1 = 1, b*x1_2*x4_3 = 1, b*x1_5*x4_2 = 1, delta*sgm*x3_0*x4_5 = 1, x1_2*x4_3 = 1, x1_6*x4_1 = 1, x1_3*x4_4 = 1, delta*sgm*x3_2*x4_1 = 1, b*x1_1*x4_7 = 1, delta*sgm*x3_4*x4_1 = 1, b*c*x1_0 = 1, delta*sgm*x3_0*x4_2 = 1, -gama*x2_1 = 1, x1_2*c = 1, -gama*sgm*x2_0*x4_7 = 1, x1_2*x4_4 = 1, b*x1_1*x4_1 = 1, -15115441682291351469493999796889194256444309189702658997423903023753904476681805540674372864157150361997078285344947881840603112801189118660874555200030360184490623136736094978269849325660690618500256904101614451942602254596539800869941050913585838514524194595778754542279355296272339717990213237690915121173658700936742021099515999046296336730983160044154731834/55966407153616880985053572623591423628062510465785284303202635187864245602512659742542807287088833506264938261188496067030530498949786507070432964268747772403978351984641552701150465895475535612767 = 1, x1_2*x4_2 = 1, x3_5 = 1, -15*gama*x4_2*sgm = 1, -alpha*x1_1 = 1, -285178703914060845731182048885212093993751879847955314363/8614621920582149333 = 1, -4*gama*sgm*x2_3*x4_1 = 1, beta*x2_1 = 1, b*x1_0*x4_7 = 1, x1_3*x4_0 = 1, x3_4 = 1, delta*sgm*x3_0*x4_0 = 1, x1_2*x4_6 = 1, x1_8*c = 1, x1_3*x4_3 = 1, x3_4*x4_1 = 1, b*c*x1_4 = 1, x1_9*x4_0 = 1, -gama = 1, beta = 1, x1_3*x4_2 = 1, -gama*sgm*x2_7*x4_0 = 1, x3_0*x4_6 = 1, beta*x2_5 = 1, -alpha*x1_2 = 1, x3_5*x4_3 = 1, delta*sgm*x3_4*x4_2 = 1, b*x1_0*x4_6 = 1, delta*sgm*x3_3*x4_0 = 1, x1_4*x4_5 = 1, b*x1_4*x4_1 = 1, -x1_8 = 1, delta*x3_1 = 1, b*x1_3*x4_0 = 1, delta*sgm*x3_2*x4_2 = 1, -gama*x2_0 = 1, -x1_9 = 1, b*c*x1_8 = 1, b*x1_0*x4_8 = 1, b*x1_6*x4_1 = 1, -5629317267317686165489021192098964322447737016015043776670822618707790786426718509586239677710877917116627651791041876428312678539839048764250576436518854083802900965418286907046686072707809546263116255395066494176408946276651216573074061161588942758327385018249058296052529508589053836325/153246013764087925164908048586910532790333032454155717482009559419018673359698329123527617565138776176515502005311383925632499363141958349866424164716281220390053 = 1, -35*gama*sgm*x2_3*x4_4 = 1, x1_3*x4_1 = 1, -4*gama*sgm*x2_1*x4_3 = 1, delta*sgm*x3_3*x4_2 = 1, -gama*x2_6 = 1, x1_1*x4_7 = 1, x1_7*x4_0 = 1, x3_2*x4_2 = 1, x3_1*x4_2 = 1, -x1_3 = 1, -10*gama*sgm*x2_3*x4_2 = 1, x3_2*x4_6 = 1, -x1_5 = 1, -gama*x2_2 = 1, delta*sgm*x3_1*x4_1 = 1, x3_3*x4_5 = 1, -10*gama*sgm*x2_2*x4_3 = 1, b*x1_2*x4_4 = 1, b*x1_4*x4_4 = 1, delta*sgm*x3_1*x4_2 = 1, -21*gama*sgm*x2_2*x4_5 = 1, x1_3*x4_6 = 1, b*x1_0*x4_0 = 1, b*x1_6*x4_2 = 1, x1_1*x4_6 = 1, x3_0*x4_3 = 1, x1_4*x4_1 = 1, -x1_0 = 1, -alpha*x1_0 = 1, delta*sgm*x3_5*x4_2 = 1, -gama*sgm*x2_1*x4_0 = 1, b*x1_1*x4_3 = 1, b*x1_1*x4_0 = 1, x3_0*x4_7 = 1, -5*gama*sgm*x2_1*x4_4 = 1, x3_4*x4_4 = 1, delta*x4_0*sgm = 1, x3_1*x4_7 = 1, b*x1_8*x4_0 = 1, b*x1_0*x4_4 = 1, x1_2*x4_5 = 1, -alpha*x1_6 = 1, x1_2*x4_7 = 1, x1_7*c = 1, b*x1_1*x4_4 = 1, x3_1*x4_3 = 1, x3_2*x4_1 = 1, -7*gama*sgm*x2_1*x4_6 = 1, b*x1_3*x4_4 = 1, x1_3*c = 1, delta*sgm*x3_2*x4_5 = 1, -gama*x2_5 = 1, b*c*x1_5 = 1, b*x1_3*x4_5 = 1, delta*sgm*x3_3*x4_3 = 1, b*x1_1*x4_5 = 1, b*x1_0*x4_3 = 1, -alpha*x1_5 = 1, -21*gama*sgm*x2_5*x4_2 = 1, b*x1_2*x4_6 = 1, x1_7*x4_1 = 1, delta = 1, -gama*sgm*x2_0*x4_5 = 1, x1_1*x4_5 = 1, b*x1_2*x4_0 = 1, b*x1_3*x4_2 = 1, -x1_1 = 1, delta*sgm*x3_1*x4_3 = 1, b*x1_3*x4_3 = 1, -6*gama*sgm*x2_2*x4_2 = 1, delta*sgm*x3_0*x4_6 = 1, b*c*x1_2 = 1, -gama*sgm*x2_5*x4_0 = 1]
# 282, -5.577999884
# 2
# [x3_6 = 6, x2_4 = 7]
# [x3_6 = [4, 2, 1, 4, 2, 2], x2_4 = [4, 1, 4, 2, 2, 4, 4]]
# [x3_6 = [1, 1, 1, 1, 1, 1], x2_4 = [1, 1, 1, 1, 1, 1, 1]]