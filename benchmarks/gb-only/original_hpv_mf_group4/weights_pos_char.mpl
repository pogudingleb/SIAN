kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
infolevel[Groebner]:=10;
et_hat:=[97508118953453359597784369735-IGM_0-IOGM_0, -IGF_0*SM_0^2*betaGGFM^2-IOF_0*SM_0^2*betaOGFM^2-IOGF_0*SM_0^2*betaGGFM^2-IOGF_0*SM_0^2*betaOGFM^2+IGF_0*IGM_0*betaGOFM^2+IGM_0*IOF_0*betaOOFM^2+IGM_0*IOGF_0*betaGOFM^2+IGM_0*IOGF_0*betaOOFM^2+IGM_0*gammaGM^2+IGM_0*mu^2+IGM_0*nuGOM^2-IOGM_0*gammaOM^2+IGM_1, -IGF_0*IGM_0*betaGOFM^2-IGF_0*IOM_0*betaGGFM^2-IGM_0*IOF_0*betaOOFM^2-IGM_0*IOGF_0*betaGOFM^2-IGM_0*IOGF_0*betaOOFM^2-IOF_0*IOM_0*betaOGFM^2-IOGF_0*IOM_0*betaGGFM^2-IOGF_0*IOM_0*betaOGFM^2-IGM_0*nuGOM^2+IOGM_0*gammaGM^2+IOGM_0*gammaOM^2+IOGM_0*mu^2-IOM_0*nuOGM^2+IOGM_1, 199245943113362595385217534663-IOM_0-IOGM_0, -IGF_0*SM_0^2*betaGOFM^2-IOF_0*SM_0^2*betaOOFM^2-IOGF_0*SM_0^2*betaGOFM^2-IOGF_0*SM_0^2*betaOOFM^2+IGF_0*IOM_0*betaGGFM^2+IOF_0*IOM_0*betaOGFM^2+IOGF_0*IOM_0*betaGGFM^2+IOGF_0*IOM_0*betaOGFM^2-IOGM_0*gammaGM^2+IOM_0*gammaOM^2+IOM_0*mu^2+IOM_0*nuOGM^2+IOM_1, 175517145482571466658437385305-IGF_0-IOGF_0, -IGM_0*SF_0^2*betaGGMF^2-IOGM_0*SF_0^2*betaGGMF^2-IOGM_0*SF_0^2*betaOGMF^2-IOM_0*SF_0^2*betaOGMF^2+IGF_0*IGM_0*betaGOMF^2+IGF_0*IOGM_0*betaGOMF^2+IGF_0*IOGM_0*betaOOMF^2+IGF_0*IOM_0*betaOOMF^2+IGF_0*gammaGF^2+IGF_0*mu^2+IGF_0*nuGOF^2-IOGF_0*gammaOF^2+IGF_1, -IGF_0*IGM_0*betaGOMF^2-IGF_0*IOGM_0*betaGOMF^2-IGF_0*IOGM_0*betaOOMF^2-IGF_0*IOM_0*betaOOMF^2-IGM_0*IOF_0*betaGGMF^2-IOF_0*IOGM_0*betaGGMF^2-IOF_0*IOGM_0*betaOGMF^2-IOF_0*IOM_0*betaOGMF^2-IGF_0*nuGOF^2-IOF_0*nuOGF^2+IOGF_0*gammaGF^2+IOGF_0*gammaOF^2+IOGF_0*mu^2+IOGF_1, 285221478830593786071161848334-IOF_0-IOGF_0, -IGM_0*SF_0^2*betaGOMF^2-IOGM_0*SF_0^2*betaGOMF^2-IOGM_0*SF_0^2*betaOOMF^2-IOM_0*SF_0^2*betaOOMF^2+IGM_0*IOF_0*betaGGMF^2+IOF_0*IOGM_0*betaGGMF^2+IOF_0*IOGM_0*betaOGMF^2+IOF_0*IOM_0*betaOGMF^2+IOF_0*gammaOF^2+IOF_0*mu^2+IOF_0*nuOGF^2-IOGF_0*gammaGF^2+IOF_1, -IGM_1-IOGM_1+7283912742264592796308753174672063618681336489577054048393694803508592524648065957588711, (nuGOM^2+gammaGM^2+mu^2+betaOOFM^2*(IOF_0+IOGF_0)+betaGOFM^2*(IGF_0+IOGF_0))*IGM_1+(-IGF_0*SM_1^2-IGF_1*SM_0^2-IOGF_0*SM_1^2-IOGF_1*SM_0^2)*betaGGFM^2+IGM_0*(IOGF_1+IGF_1)*betaGOFM^2+(-IOF_0*SM_1^2-IOF_1*SM_0^2-IOGF_0*SM_1^2-IOGF_1*SM_0^2)*betaOGFM^2+IGM_0*(IOGF_1+IOF_1)*betaOOFM^2-gammaOM^2*IOGM_1+IGM_2, ((-IGF_0-IOGF_0)*betaGOFM^2+(-IOF_0-IOGF_0)*betaOOFM^2-nuGOM^2)*IGM_1+((-IGF_0-IOGF_0)*betaGGFM^2+(-IOF_0-IOGF_0)*betaOGFM^2-nuOGM^2)*IOM_1-IOM_0*(IOGF_1+IGF_1)*betaGGFM^2-IGM_0*(IOGF_1+IGF_1)*betaGOFM^2-IOM_0*(IOGF_1+IOF_1)*betaOGFM^2-IGM_0*(IOGF_1+IOF_1)*betaOOFM^2+(gammaGM^2+gammaOM^2+mu^2)*IOGM_1+IOGM_2, 2*IGF_0*SM_0^2*betaGGFM^2+2*IGF_0*SM_0^2*betaGOFM^2+2*IOF_0*SM_0^2*betaOGFM^2+2*IOF_0*SM_0^2*betaOOFM^2+2*IOGF_0*SM_0^2*betaGGFM^2+2*IOGF_0*SM_0^2*betaGOFM^2+2*IOGF_0*SM_0^2*betaOGFM^2+2*IOGF_0*SM_0^2*betaOOFM^2+2*SM_0^2*mu^2-2*IGM_0*gammaGM^2-2*IOM_0*gammaOM^2+2*SM_1^2-mu^2, -IOGM_1-IOM_1+1623528416966471363773867264441616509415163364884026520175525654435376507872121018167651, (nuOGM^2+gammaOM^2+mu^2+betaOGFM^2*(IOF_0+IOGF_0)+betaGGFM^2*(IGF_0+IOGF_0))*IOM_1+IOM_0*(IOGF_1+IGF_1)*betaGGFM^2+(-IGF_0*SM_1^2-IGF_1*SM_0^2-IOGF_0*SM_1^2-IOGF_1*SM_0^2)*betaGOFM^2+IOM_0*(IOGF_1+IOF_1)*betaOGFM^2+(-IOF_0*SM_1^2-IOF_1*SM_0^2-IOGF_0*SM_1^2-IOGF_1*SM_0^2)*betaOOFM^2-gammaGM^2*IOGM_1+IOM_2, -IGF_1-IOGF_1+3693948057289645380986651123185169114270313216174342100113020784323093172112839890455451, (nuGOF^2+gammaGF^2+mu^2+betaOOMF^2*(IOM_0+IOGM_0)+betaGOMF^2*(IGM_0+IOGM_0))*IGF_1+(-IGM_0*SF_1^2-IGM_1*SF_0^2-IOGM_0*SF_1^2-IOGM_1*SF_0^2)*betaGGMF^2+IGF_0*(IOGM_1+IGM_1)*betaGOMF^2+(-IOGM_0*SF_1^2-IOGM_1*SF_0^2-IOM_0*SF_1^2-IOM_1*SF_0^2)*betaOGMF^2+IGF_0*(IOM_1+IOGM_1)*betaOOMF^2-gammaOF^2*IOGF_1+IGF_2, ((-IGM_0-IOGM_0)*betaGOMF^2+(-IOGM_0-IOM_0)*betaOOMF^2-nuGOF^2)*IGF_1+((-IGM_0-IOGM_0)*betaGGMF^2+(-IOGM_0-IOM_0)*betaOGMF^2-nuOGF^2)*IOF_1-IOF_0*(IOGM_1+IGM_1)*betaGGMF^2-IGF_0*(IOGM_1+IGM_1)*betaGOMF^2-IOF_0*(IOM_1+IOGM_1)*betaOGMF^2-IGF_0*(IOM_1+IOGM_1)*betaOOMF^2+(gammaGF^2+gammaOF^2+mu^2)*IOGF_1+IOGF_2, 2*IGM_0*SF_0^2*betaGGMF^2+2*IGM_0*SF_0^2*betaGOMF^2+2*IOGM_0*SF_0^2*betaGGMF^2+2*IOGM_0*SF_0^2*betaGOMF^2+2*IOGM_0*SF_0^2*betaOGMF^2+2*IOGM_0*SF_0^2*betaOOMF^2+2*IOM_0*SF_0^2*betaOGMF^2+2*IOM_0*SF_0^2*betaOOMF^2+2*SF_0^2*mu^2-2*IGF_0*gammaGF^2-2*IOF_0*gammaOF^2+2*SF_1^2-mu^2, -IOF_1-IOGF_1+1422466831033510637161257209097378483511007890241624970508560220331737129492657922208634, (nuOGF^2+gammaOF^2+mu^2+betaOGMF^2*(IOM_0+IOGM_0)+betaGGMF^2*(IGM_0+IOGM_0))*IOF_1+IOF_0*(IOGM_1+IGM_1)*betaGGMF^2+(-IGM_0*SF_1^2-IGM_1*SF_0^2-IOGM_0*SF_1^2-IOGM_1*SF_0^2)*betaGOMF^2+IOF_0*(IOM_1+IOGM_1)*betaOGMF^2+(-IOGM_0*SF_1^2-IOGM_1*SF_0^2-IOM_0*SF_1^2-IOM_1*SF_0^2)*betaOOMF^2-gammaGF^2*IOGF_1+IOF_2, -IGM_2-IOGM_2-327131288738826187785649302776524543450421934574842240234930231801426072058372615695822479309110667749984242025823252268653780197734893332620142276, (nuGOM^2+gammaGM^2+mu^2+betaOOFM^2*(IOF_0+IOGF_0)+betaGOFM^2*(IGF_0+IOGF_0))*IGM_2+(-IGF_0*SM_2^2-2*IGF_1*SM_1^2-IGF_2*SM_0^2-IOGF_0*SM_2^2-2*IOGF_1*SM_1^2-IOGF_2*SM_0^2)*betaGGFM^2+(2*IGF_1*IGM_1+IGF_2*IGM_0+IGM_0*IOGF_2+2*IGM_1*IOGF_1)*betaGOFM^2+(-IOF_0*SM_2^2-2*IOF_1*SM_1^2-IOF_2*SM_0^2-IOGF_0*SM_2^2-2*IOGF_1*SM_1^2-IOGF_2*SM_0^2)*betaOGFM^2+(IGM_0*IOF_2+IGM_0*IOGF_2+2*IGM_1*IOF_1+2*IGM_1*IOGF_1)*betaOOFM^2-gammaOM^2*IOGM_2+IGM_3, ((-IGF_0-IOGF_0)*IOM_2-IOM_0*IOGF_2-IOM_0*IGF_2-2*IOM_1*IOGF_1-2*IOM_1*IGF_1)*betaGGFM^2+((-IGF_0-IOGF_0)*IGM_2-2*IOGF_1*IGM_1-IOGF_2*IGM_0-2*IGF_1*IGM_1-IGF_2*IGM_0)*betaGOFM^2+((-IOF_0-IOGF_0)*IOM_2-IOF_2*IOM_0-IOM_0*IOGF_2-2*IOM_1*IOGF_1-2*IOM_1*IOF_1)*betaOGFM^2+((-IOF_0-IOGF_0)*IGM_2-IOF_2*IGM_0-2*IOGF_1*IGM_1-IOGF_2*IGM_0-2*IGM_1*IOF_1)*betaOOFM^2-IGM_2*nuGOM^2-IOM_2*nuOGM^2+(gammaGM^2+gammaOM^2+mu^2)*IOGM_2+IOGM_3, ((2*IOGF_0+2*IGF_0)*betaGGFM^2+(2*IOGF_0+2*IGF_0)*betaGOFM^2+(2*IOGF_0+2*IOF_0)*betaOGFM^2+(2*IOGF_0+2*IOF_0)*betaOOFM^2+2*mu^2)*SM_1^2+((2*IOGF_1+2*IGF_1)*betaGGFM^2+(2*IOGF_1+2*IGF_1)*betaGOFM^2+2*(betaOGFM^2+betaOOFM^2)*(IOGF_1+IOF_1))*SM_0^2-2*gammaOM^2*IOM_1-2*gammaGM^2*IGM_1+2*SM_2^2, -IOGM_2-IOM_2-123213868372558311470475380955602985797873726388615018586799171367203963929927483609900907307852534003794475283423023203358111267589978160598104369/2, (nuOGM^2+gammaOM^2+mu^2+betaOGFM^2*(IOF_0+IOGF_0)+betaGGFM^2*(IGF_0+IOGF_0))*IOM_2+(2*IGF_1*IOM_1+IGF_2*IOM_0+2*IOGF_1*IOM_1+IOGF_2*IOM_0)*betaGGFM^2+(-IGF_0*SM_2^2-2*IGF_1*SM_1^2-IGF_2*SM_0^2-IOGF_0*SM_2^2-2*IOGF_1*SM_1^2-IOGF_2*SM_0^2)*betaGOFM^2+(2*IOF_1*IOM_1+IOF_2*IOM_0+2*IOGF_1*IOM_1+IOGF_2*IOM_0)*betaOGFM^2+(-IOF_0*SM_2^2-2*IOF_1*SM_1^2-IOF_2*SM_0^2-IOGF_0*SM_2^2-2*IOGF_1*SM_1^2-IOGF_2*SM_0^2)*betaOOFM^2-gammaGM^2*IOGM_2+IOM_3, -IGF_2-IOGF_2+295373479984315807925951905952530319057929015347472895106640867838990528880623717516678138117429809921193697604422861525471912391653102447136957369/2, (nuGOF^2+gammaGF^2+mu^2+betaOOMF^2*(IOM_0+IOGM_0)+betaGOMF^2*(IGM_0+IOGM_0))*IGF_2+(-IGM_0*SF_2^2-2*IGM_1*SF_1^2-IGM_2*SF_0^2-IOGM_0*SF_2^2-2*IOGM_1*SF_1^2-IOGM_2*SF_0^2)*betaGGMF^2+(IGF_0*IGM_2+IGF_0*IOGM_2+2*IGF_1*IGM_1+2*IGF_1*IOGM_1)*betaGOMF^2+(-IOGM_0*SF_2^2-2*IOGM_1*SF_1^2-IOGM_2*SF_0^2-IOM_0*SF_2^2-2*IOM_1*SF_1^2-IOM_2*SF_0^2)*betaOGMF^2+(IGF_0*IOGM_2+IGF_0*IOM_2+2*IGF_1*IOGM_1+2*IGF_1*IOM_1)*betaOOMF^2-gammaOF^2*IOGF_2+IGF_3, ((-IGM_0-IOGM_0)*IOF_2-2*IOGM_1*IOF_1-IOGM_2*IOF_0-2*IGM_1*IOF_1-IGM_2*IOF_0)*betaGGMF^2+((-IGM_0-IOGM_0)*IGF_2-2*IOGM_1*IGF_1-IGF_0*IOGM_2-IGF_0*IGM_2-2*IGF_1*IGM_1)*betaGOMF^2+((-IOGM_0-IOM_0)*IOF_2-2*IOM_1*IOF_1-IOM_2*IOF_0-2*IOGM_1*IOF_1-IOGM_2*IOF_0)*betaOGMF^2+((-IOGM_0-IOM_0)*IGF_2-2*IOM_1*IGF_1-IOM_2*IGF_0-2*IOGM_1*IGF_1-IGF_0*IOGM_2)*betaOOMF^2-IGF_2*nuGOF^2-IOF_2*nuOGF^2+(gammaGF^2+gammaOF^2+mu^2)*IOGF_2+IOGF_3, ((2*IOGM_0+2*IGM_0)*betaGGMF^2+(2*IOGM_0+2*IGM_0)*betaGOMF^2+(2*IOGM_0+2*IOM_0)*betaOGMF^2+(2*IOGM_0+2*IOM_0)*betaOOMF^2+2*mu^2)*SF_1^2+((2*IOGM_1+2*IGM_1)*betaGGMF^2+(2*IOGM_1+2*IGM_1)*betaGOMF^2+2*(betaOGMF^2+betaOOMF^2)*(IOM_1+IOGM_1))*SF_0^2-2*gammaGF^2*IGF_1-2*gammaOF^2*IOF_1+2*SF_2^2, -IOF_2-IOGF_2+107776252533508028276946369092821768482763737435504460465883228048399079625745952349941748451013544387104147434323587413275958028884257557885332911/2, (nuOGF^2+gammaOF^2+mu^2+betaOGMF^2*(IOM_0+IOGM_0)+betaGGMF^2*(IGM_0+IOGM_0))*IOF_2+(2*IGM_1*IOF_1+IGM_2*IOF_0+IOF_0*IOGM_2+2*IOF_1*IOGM_1)*betaGGMF^2+(-IGM_0*SF_2^2-2*IGM_1*SF_1^2-IGM_2*SF_0^2-IOGM_0*SF_2^2-2*IOGM_1*SF_1^2-IOGM_2*SF_0^2)*betaGOMF^2+(IOF_0*IOGM_2+IOF_0*IOM_2+2*IOF_1*IOGM_1+2*IOF_1*IOM_1)*betaOGMF^2+(-IOGM_0*SF_2^2-2*IOGM_1*SF_1^2-IOGM_2*SF_0^2-IOM_0*SF_2^2-2*IOM_1*SF_1^2-IOM_2*SF_0^2)*betaOOMF^2-gammaGF^2*IOGF_2+IOF_3, -IGM_3-IOGM_3+26514236107588614635621237751190876307659947432862864433849198508090320599221113379888826729209093953826928520081848027393803226170825817320145292066625466077889994278003743538156935065678567132835131541733/2, (-IGF_0*SM_3^2-3*IGF_1*SM_2^2-3*IGF_2*SM_1^2-IGF_3*SM_0^2-IOGF_0*SM_3^2-3*IOGF_1*SM_2^2-3*IOGF_2*SM_1^2-IOGF_3*SM_0^2)*betaGGFM^2+((IGF_0+IOGF_0)*IGM_3+IGM_0*IOGF_3+IGM_0*IGF_3+3*IGM_1*IOGF_2+3*IGM_1*IGF_2+3*IGM_2*IOGF_1+3*IGM_2*IGF_1)*betaGOFM^2+(-IOF_0*SM_3^2-3*IOF_1*SM_2^2-3*IOF_2*SM_1^2-IOF_3*SM_0^2-IOGF_0*SM_3^2-3*IOGF_1*SM_2^2-3*IOGF_2*SM_1^2-IOGF_3*SM_0^2)*betaOGFM^2+((IOF_0+IOGF_0)*IGM_3+IGM_0*IOF_3+IGM_0*IOGF_3+3*IGM_1*IOF_2+3*IGM_1*IOGF_2+3*IGM_2*IOF_1+3*IGM_2*IOGF_1)*betaOOFM^2+(gammaGM^2+mu^2+nuGOM^2)*IGM_3-gammaOM^2*IOGM_3+IGM_4, ((-IGF_0-IOGF_0)*IOM_3-IOM_0*IOGF_3-IOM_0*IGF_3-3*IOM_1*IOGF_2-3*IOM_1*IGF_2-3*IOM_2*IOGF_1-3*IOM_2*IGF_1)*betaGGFM^2+((-IGF_0-IOGF_0)*IGM_3-3*IGM_2*IOGF_1-3*IGM_1*IOGF_2-IGM_0*IOGF_3-3*IGM_2*IGF_1-3*IGM_1*IGF_2-IGM_0*IGF_3)*betaGOFM^2+((-IOF_0-IOGF_0)*IOM_3-3*IOM_2*IOF_1-3*IOF_2*IOM_1-IOF_3*IOM_0-IOM_0*IOGF_3-3*IOM_1*IOGF_2-3*IOM_2*IOGF_1)*betaOGFM^2+((-IOF_0-IOGF_0)*IGM_3-3*IGM_2*IOF_1-3*IGM_1*IOF_2-IGM_0*IOF_3-3*IGM_2*IOGF_1-3*IGM_1*IOGF_2-IGM_0*IOGF_3)*betaOOFM^2-IGM_3*nuGOM^2-IOM_3*nuOGM^2+(gammaGM^2+gammaOM^2+mu^2)*IOGM_3+IOGM_4, ((2*IOGF_0+2*IGF_0)*betaGGFM^2+(2*IOGF_0+2*IGF_0)*betaGOFM^2+(2*IOGF_0+2*IOF_0)*betaOGFM^2+(2*IOGF_0+2*IOF_0)*betaOOFM^2+2*mu^2)*SM_2^2+((2*IOGF_2+2*IGF_2)*betaGGFM^2+(2*IOGF_2+2*IGF_2)*betaGOFM^2+2*(betaOGFM^2+betaOOFM^2)*(IOF_2+IOGF_2))*SM_0^2+((4*IOGF_1+4*IGF_1)*betaGGFM^2+(4*IOGF_1+4*IGF_1)*betaGOFM^2+4*(betaOGFM^2+betaOOFM^2)*(IOGF_1+IOF_1))*SM_1^2-2*gammaOM^2*IOM_2-2*gammaGM^2*IGM_2+2*SM_3^2, -IOGM_3-IOM_3+2083817365208573936255208542674340724326149601061590154110453228655623823799130230700773323431261958466897976938072757051504005881264385980933903536758068901150949845510132750946982332864042073521435951334, ((IGF_0+IOGF_0)*IOM_3+IOM_0*IOGF_3+IOM_0*IGF_3+3*IOM_1*IOGF_2+3*IOM_1*IGF_2+3*IOM_2*IOGF_1+3*IOM_2*IGF_1)*betaGGFM^2+(-IGF_0*SM_3^2-3*IGF_1*SM_2^2-3*IGF_2*SM_1^2-IGF_3*SM_0^2-IOGF_0*SM_3^2-3*IOGF_1*SM_2^2-3*IOGF_2*SM_1^2-IOGF_3*SM_0^2)*betaGOFM^2+((IOF_0+IOGF_0)*IOM_3+3*IOM_2*IOF_1+3*IOF_2*IOM_1+IOF_3*IOM_0+IOM_0*IOGF_3+3*IOM_1*IOGF_2+3*IOM_2*IOGF_1)*betaOGFM^2+(-IOF_0*SM_3^2-3*IOF_1*SM_2^2-3*IOF_2*SM_1^2-IOF_3*SM_0^2-IOGF_0*SM_3^2-3*IOGF_1*SM_2^2-3*IOGF_2*SM_1^2-IOGF_3*SM_0^2)*betaOOFM^2+(gammaOM^2+mu^2+nuOGM^2)*IOM_3-gammaGM^2*IOGM_3+IOM_4, -IGF_3-IOGF_3-41266469183115959522799806046415427562030781253152399748324977525975161448143799255557694961514744071584442662829404200695792407996155070723540803388756905046793633227426399246761816213087100914868887990603/2, (-IGM_0*SF_3^2-3*IGM_1*SF_2^2-3*IGM_2*SF_1^2-IGM_3*SF_0^2-IOGM_0*SF_3^2-3*IOGM_1*SF_2^2-3*IOGM_2*SF_1^2-IOGM_3*SF_0^2)*betaGGMF^2+((IGM_0+IOGM_0)*IGF_3+3*IGM_1*IGF_2+3*IGM_2*IGF_1+IGM_3*IGF_0+3*IOGM_1*IGF_2+3*IOGM_2*IGF_1+IOGM_3*IGF_0)*betaGOMF^2+(-IOGM_0*SF_3^2-3*IOGM_1*SF_2^2-3*IOGM_2*SF_1^2-IOGM_3*SF_0^2-IOM_0*SF_3^2-3*IOM_1*SF_2^2-3*IOM_2*SF_1^2-IOM_3*SF_0^2)*betaOGMF^2+((IOM_0+IOGM_0)*IGF_3+3*IOM_1*IGF_2+3*IOM_2*IGF_1+IOM_3*IGF_0+3*IOGM_1*IGF_2+3*IOGM_2*IGF_1+IOGM_3*IGF_0)*betaOOMF^2+(gammaGF^2+mu^2+nuGOF^2)*IGF_3-gammaOF^2*IOGF_3+IGF_4, ((-IGM_0-IOGM_0)*IOF_3-IOF_0*IOGM_3-IOF_0*IGM_3-3*IOF_1*IOGM_2-3*IGM_2*IOF_1-3*IOF_2*IOGM_1-3*IGM_1*IOF_2)*betaGGMF^2+((-IGM_0-IOGM_0)*IGF_3-IOGM_3*IGF_0-IGM_3*IGF_0-3*IOGM_1*IGF_2-3*IOGM_2*IGF_1-3*IGM_2*IGF_1-3*IGM_1*IGF_2)*betaGOMF^2+((-IOGM_0-IOM_0)*IOF_3-IOM_3*IOF_0-IOF_0*IOGM_3-3*IOM_2*IOF_1-3*IOF_1*IOGM_2-3*IOF_2*IOM_1-3*IOF_2*IOGM_1)*betaOGMF^2+((-IOGM_0-IOM_0)*IGF_3-3*IOM_1*IGF_2-3*IOM_2*IGF_1-IOM_3*IGF_0-IOGM_3*IGF_0-3*IOGM_1*IGF_2-3*IOGM_2*IGF_1)*betaOOMF^2-IGF_3*nuGOF^2-IOF_3*nuOGF^2+(gammaGF^2+gammaOF^2+mu^2)*IOGF_3+IOGF_4, ((2*IOGM_0+2*IGM_0)*betaGGMF^2+(2*IOGM_0+2*IGM_0)*betaGOMF^2+(2*IOGM_0+2*IOM_0)*betaOGMF^2+(2*IOGM_0+2*IOM_0)*betaOOMF^2+2*mu^2)*SF_2^2+((2*IOGM_2+2*IGM_2)*betaGGMF^2+(2*IOGM_2+2*IGM_2)*betaGOMF^2+2*(betaOGMF^2+betaOOMF^2)*(IOM_2+IOGM_2))*SF_0^2+((4*IOGM_1+4*IGM_1)*betaGGMF^2+(4*IOGM_1+4*IGM_1)*betaGOMF^2+4*(betaOGMF^2+betaOOMF^2)*(IOM_1+IOGM_1))*SF_1^2-2*gammaOF^2*IOF_2-2*gammaGF^2*IGF_2+2*SF_3^2, -IOF_3-IOGF_3-6850354874195166649707588809223335872617714261928598554541276598639489331192277137276139613890633092746573427646651925707926127344668500267284396748785696145746623752184236456807776554494914436135467941396, ((IGM_0+IOGM_0)*IOF_3+IOF_0*IOGM_3+IOF_0*IGM_3+3*IOF_1*IOGM_2+3*IGM_2*IOF_1+3*IOF_2*IOGM_1+3*IGM_1*IOF_2)*betaGGMF^2+(-IGM_0*SF_3^2-3*IGM_1*SF_2^2-3*IGM_2*SF_1^2-IGM_3*SF_0^2-IOGM_0*SF_3^2-3*IOGM_1*SF_2^2-3*IOGM_2*SF_1^2-IOGM_3*SF_0^2)*betaGOMF^2+((IOM_0+IOGM_0)*IOF_3+IOM_3*IOF_0+IOF_0*IOGM_3+3*IOM_2*IOF_1+3*IOF_1*IOGM_2+3*IOF_2*IOM_1+3*IOF_2*IOGM_1)*betaOGMF^2+(-IOGM_0*SF_3^2-3*IOGM_1*SF_2^2-3*IOGM_2*SF_1^2-IOGM_3*SF_0^2-IOM_0*SF_3^2-3*IOM_1*SF_2^2-3*IOM_2*SF_1^2-IOM_3*SF_0^2)*betaOOMF^2+(gammaOF^2+mu^2+nuOGF^2)*IOF_3-gammaGF^2*IOGF_3+IOF_4, -IGM_4-IOGM_4-1953424670071179430591412762051547840076116509938515302727168138176932383790970521917056172242798291420534952707278914115190844783629599398487260909296018014085274948714764224571206216726276247549299776998774905579199342423243664667616142561472640040118865008696655/2, (-IGF_0*SM_4^2-4*IGF_1*SM_3^2-6*IGF_2*SM_2^2-4*IGF_3*SM_1^2-IGF_4*SM_0^2-IOGF_0*SM_4^2-4*IOGF_1*SM_3^2-6*IOGF_2*SM_2^2-4*IOGF_3*SM_1^2-IOGF_4*SM_0^2)*betaGGFM^2+((IGF_0+IOGF_0)*IGM_4+4*IGF_3*IGM_1+IGF_4*IGM_0+IGM_0*IOGF_4+4*IGM_1*IOGF_3+6*IGM_2*IOGF_2+6*IGF_2*IGM_2+4*IGM_3*IOGF_1+4*IGF_1*IGM_3)*betaGOFM^2+(-IOF_0*SM_4^2-4*IOF_1*SM_3^2-6*IOF_2*SM_2^2-4*IOF_3*SM_1^2-IOF_4*SM_0^2-IOGF_0*SM_4^2-4*IOGF_1*SM_3^2-6*IOGF_2*SM_2^2-4*IOGF_3*SM_1^2-IOGF_4*SM_0^2)*betaOGFM^2+((IOF_0+IOGF_0)*IGM_4+IGM_0*IOF_4+IGM_0*IOGF_4+4*IGM_1*IOF_3+4*IGM_1*IOGF_3+6*IGM_2*IOF_2+6*IGM_2*IOGF_2+4*IGM_3*IOF_1+4*IGM_3*IOGF_1)*betaOOFM^2+(gammaGM^2+mu^2+nuGOM^2)*IGM_4-gammaOM^2*IOGM_4+IGM_5, ((-IGF_0-IOGF_0)*IOM_4-IOGF_4*IOM_0-IOM_0*IGF_4-4*IOGF_3*IOM_1-4*IOM_1*IGF_3-6*IOGF_2*IOM_2-6*IOM_2*IGF_2-4*IOGF_1*IOM_3-4*IOM_3*IGF_1)*betaGGFM^2+((-IGF_0-IOGF_0)*IGM_4-4*IGM_3*IOGF_1-4*IGF_1*IGM_3-6*IGM_2*IOGF_2-4*IGM_1*IOGF_3-IGM_0*IOGF_4-6*IGF_2*IGM_2-4*IGF_3*IGM_1-IGF_4*IGM_0)*betaGOFM^2+((-IOF_0-IOGF_0)*IOM_4-4*IOF_1*IOM_3-6*IOF_2*IOM_2-4*IOF_3*IOM_1-IOF_4*IOM_0-IOGF_4*IOM_0-4*IOGF_3*IOM_1-6*IOGF_2*IOM_2-4*IOGF_1*IOM_3)*betaOGFM^2+((-IOF_0-IOGF_0)*IGM_4-4*IGM_3*IOF_1-4*IGM_3*IOGF_1-6*IGM_2*IOF_2-4*IGM_1*IOF_3-IGM_0*IOF_4-6*IGM_2*IOGF_2-4*IGM_1*IOGF_3-IGM_0*IOGF_4)*betaOOFM^2-IGM_4*nuGOM^2-IOM_4*nuOGM^2+(gammaGM^2+gammaOM^2+mu^2)*IOGM_4+IOGM_5, ((2*IOGF_0+2*IGF_0)*betaGGFM^2+(2*IOGF_0+2*IGF_0)*betaGOFM^2+(2*IOGF_0+2*IOF_0)*betaOGFM^2+(2*IOGF_0+2*IOF_0)*betaOOFM^2+2*mu^2)*SM_3^2+((2*IGF_3+2*IOGF_3)*SM_0^2+(6*IOGF_2+6*IGF_2)*SM_1^2+6*SM_2^2*(IOGF_1+IGF_1))*betaGGFM^2+((2*IGF_3+2*IOGF_3)*SM_0^2+(6*IOGF_2+6*IGF_2)*SM_1^2+6*SM_2^2*(IOGF_1+IGF_1))*betaGOFM^2+((2*IOF_3+2*IOGF_3)*SM_0^2+(6*IOF_2+6*IOGF_2)*SM_1^2+6*SM_2^2*(IOGF_1+IOF_1))*betaOGFM^2+((2*IOF_3+2*IOGF_3)*SM_0^2+(6*IOF_2+6*IOGF_2)*SM_1^2+6*SM_2^2*(IOGF_1+IOF_1))*betaOOFM^2-2*gammaGM^2*IGM_3-2*gammaOM^2*IOM_3+2*SM_4^2, -IOGM_4-IOM_4-715729595451857020377097255002195333082991210294244204198241350436185786857803472977097759452259536936089077540560862078147662383687167653443208394262173449222637346293295528106656859855416537367407385316566803807316646757194876754057710157636014108324498588871837/4, ((IGF_0+IOGF_0)*IOM_4+4*IOM_1*IGF_3+IOM_0*IGF_4+IOGF_4*IOM_0+4*IOGF_3*IOM_1+6*IOGF_2*IOM_2+6*IOM_2*IGF_2+4*IOGF_1*IOM_3+4*IOM_3*IGF_1)*betaGGFM^2+(-IGF_0*SM_4^2-4*IGF_1*SM_3^2-6*IGF_2*SM_2^2-4*IGF_3*SM_1^2-IGF_4*SM_0^2-IOGF_0*SM_4^2-4*IOGF_1*SM_3^2-6*IOGF_2*SM_2^2-4*IOGF_3*SM_1^2-IOGF_4*SM_0^2)*betaGOFM^2+((IOF_0+IOGF_0)*IOM_4+4*IOF_1*IOM_3+6*IOF_2*IOM_2+4*IOF_3*IOM_1+IOF_4*IOM_0+IOGF_4*IOM_0+4*IOGF_3*IOM_1+6*IOGF_2*IOM_2+4*IOGF_1*IOM_3)*betaOGFM^2+(-IOF_0*SM_4^2-4*IOF_1*SM_3^2-6*IOF_2*SM_2^2-4*IOF_3*SM_1^2-IOF_4*SM_0^2-IOGF_0*SM_4^2-4*IOGF_1*SM_3^2-6*IOGF_2*SM_2^2-4*IOGF_3*SM_1^2-IOGF_4*SM_0^2)*betaOOFM^2+(gammaOM^2+mu^2+nuOGM^2)*IOM_4-gammaGM^2*IOGM_4+IOM_5, -IGF_4-IOGF_4+3345860007716618827512318586110565794976523787023023055711271870122105098208918471227016743159380163470935130490150804319173643025266516466415278562467049404388019434896602412659468436729359178689144662329328790683260253770057681207091611422278763847559764384276091/4, (-IGM_0*SF_4^2-4*IGM_1*SF_3^2-6*IGM_2*SF_2^2-4*IGM_3*SF_1^2-IGM_4*SF_0^2-IOGM_0*SF_4^2-4*IOGM_1*SF_3^2-6*IOGM_2*SF_2^2-4*IOGM_3*SF_1^2-IOGM_4*SF_0^2)*betaGGMF^2+((IGM_0+IOGM_0)*IGF_4+4*IGF_3*IGM_1+4*IGF_3*IOGM_1+6*IGF_2*IGM_2+4*IGF_1*IGM_3+IGM_4*IGF_0+6*IGF_2*IOGM_2+4*IGF_1*IOGM_3+IGF_0*IOGM_4)*betaGOMF^2+(-IOGM_0*SF_4^2-4*IOGM_1*SF_3^2-6*IOGM_2*SF_2^2-4*IOGM_3*SF_1^2-IOGM_4*SF_0^2-IOM_0*SF_4^2-4*IOM_1*SF_3^2-6*IOM_2*SF_2^2-4*IOM_3*SF_1^2-IOM_4*SF_0^2)*betaOGMF^2+((IOM_0+IOGM_0)*IGF_4+4*IOM_1*IGF_3+4*IGF_3*IOGM_1+6*IOM_2*IGF_2+4*IOM_3*IGF_1+IOM_4*IGF_0+6*IGF_2*IOGM_2+4*IGF_1*IOGM_3+IGF_0*IOGM_4)*betaOOMF^2+(gammaGF^2+mu^2+nuGOF^2)*IGF_4-gammaOF^2*IOGF_4+IGF_5, ((-IGM_0-IOGM_0)*IOF_4-4*IGM_3*IOF_1-IGM_4*IOF_0-IOF_0*IOGM_4-4*IOF_1*IOGM_3-6*IOF_2*IOGM_2-6*IGM_2*IOF_2-4*IOF_3*IOGM_1-4*IGM_1*IOF_3)*betaGGMF^2+((-IGM_0-IOGM_0)*IGF_4-4*IGF_1*IGM_3-IGM_4*IGF_0-IGF_0*IOGM_4-4*IGF_3*IOGM_1-6*IGF_2*IOGM_2-4*IGF_1*IOGM_3-6*IGF_2*IGM_2-4*IGF_3*IGM_1)*betaGOMF^2+((-IOGM_0-IOM_0)*IOF_4-IOF_0*IOM_4-IOF_0*IOGM_4-4*IOF_1*IOM_3-4*IOF_1*IOGM_3-6*IOF_2*IOM_2-6*IOF_2*IOGM_2-4*IOF_3*IOM_1-4*IOF_3*IOGM_1)*betaOGMF^2+((-IOGM_0-IOM_0)*IGF_4-4*IOM_1*IGF_3-6*IOM_2*IGF_2-IOM_4*IGF_0-IGF_0*IOGM_4-4*IOM_3*IGF_1-4*IGF_3*IOGM_1-6*IGF_2*IOGM_2-4*IGF_1*IOGM_3)*betaOOMF^2-IGF_4*nuGOF^2-IOF_4*nuOGF^2+(gammaGF^2+gammaOF^2+mu^2)*IOGF_4+IOGF_5, ((2*IOGM_0+2*IGM_0)*betaGGMF^2+(2*IOGM_0+2*IGM_0)*betaGOMF^2+(2*IOGM_0+2*IOM_0)*betaOGMF^2+(2*IOGM_0+2*IOM_0)*betaOOMF^2+2*mu^2)*SF_3^2+((2*IGM_3+2*IOGM_3)*SF_0^2+(6*IGM_2+6*IOGM_2)*SF_1^2+6*SF_2^2*(IOGM_1+IGM_1))*betaGGMF^2+((2*IGM_3+2*IOGM_3)*SF_0^2+(6*IGM_2+6*IOGM_2)*SF_1^2+6*SF_2^2*(IOGM_1+IGM_1))*betaGOMF^2+((2*IOM_3+2*IOGM_3)*SF_0^2+(6*IOM_2+6*IOGM_2)*SF_1^2+6*SF_2^2*(IOM_1+IOGM_1))*betaOGMF^2+((2*IOM_3+2*IOGM_3)*SF_0^2+(6*IOM_2+6*IOGM_2)*SF_1^2+6*SF_2^2*(IOM_1+IOGM_1))*betaOOMF^2-2*gammaOF^2*IOF_3-2*gammaGF^2*IGF_3+2*SF_4^2, -IOF_4-IOGF_4+552241057006457457905375098660764124284592299014019977463409361312054016906132656537385170231622307783895227121958291986897530632167586743388967749488829180269269985226869191534450901909129377639364095809481386672632247440491076838300102732486040279501471112206013/2, ((IGM_0+IOGM_0)*IOF_4+4*IGM_3*IOF_1+IGM_4*IOF_0+IOF_0*IOGM_4+4*IOF_1*IOGM_3+6*IOF_2*IOGM_2+6*IGM_2*IOF_2+4*IOF_3*IOGM_1+4*IGM_1*IOF_3)*betaGGMF^2+(-IGM_0*SF_4^2-4*IGM_1*SF_3^2-6*IGM_2*SF_2^2-4*IGM_3*SF_1^2-IGM_4*SF_0^2-IOGM_0*SF_4^2-4*IOGM_1*SF_3^2-6*IOGM_2*SF_2^2-4*IOGM_3*SF_1^2-IOGM_4*SF_0^2)*betaGOMF^2+((IOM_0+IOGM_0)*IOF_4+IOF_0*IOM_4+IOF_0*IOGM_4+4*IOF_1*IOM_3+4*IOF_1*IOGM_3+6*IOF_2*IOM_2+6*IOF_2*IOGM_2+4*IOF_3*IOM_1+4*IOF_3*IOGM_1)*betaOGMF^2+(-IOGM_0*SF_4^2-4*IOGM_1*SF_3^2-6*IOGM_2*SF_2^2-4*IOGM_3*SF_1^2-IOGM_4*SF_0^2-IOM_0*SF_4^2-4*IOM_1*SF_3^2-6*IOM_2*SF_2^2-4*IOM_3*SF_1^2-IOM_4*SF_0^2)*betaOOMF^2+(gammaOF^2+mu^2+nuOGF^2)*IOF_4-gammaGF^2*IOGF_4+IOF_5, -IGM_5-IOGM_5+584409042182716010092405599813116126497311644205724682647722994644372308434958266281899270190677948311504964192034648598359221660415470670580487900619578538017213257310353481210278493425501113696252572592991890525392002160246044396689550590917662386044592817587011969791349008306771667360905718711744834356201097373938722919/4, (-IGF_0*SM_5^2-5*IGF_1*SM_4^2-10*IGF_2*SM_3^2-10*IGF_3*SM_2^2-5*IGF_4*SM_1^2-IGF_5*SM_0^2-IOGF_0*SM_5^2-5*IOGF_1*SM_4^2-10*IOGF_2*SM_3^2-10*IOGF_3*SM_2^2-5*IOGF_4*SM_1^2-IOGF_5*SM_0^2)*betaGGFM^2+((IGF_0+IOGF_0)*IGM_5+10*IGF_2*IGM_3+10*IGF_3*IGM_2+5*IGF_4*IGM_1+IGF_5*IGM_0+IGM_0*IOGF_5+5*IGM_1*IOGF_4+10*IGM_2*IOGF_3+10*IGM_3*IOGF_2+5*IGM_4*IOGF_1+5*IGF_1*IGM_4)*betaGOFM^2+(-IOF_0*SM_5^2-5*IOF_1*SM_4^2-10*IOF_2*SM_3^2-10*IOF_3*SM_2^2-5*IOF_4*SM_1^2-IOF_5*SM_0^2-IOGF_0*SM_5^2-5*IOGF_1*SM_4^2-10*IOGF_2*SM_3^2-10*IOGF_3*SM_2^2-5*IOGF_4*SM_1^2-IOGF_5*SM_0^2)*betaOGFM^2+((IOF_0+IOGF_0)*IGM_5+IGM_0*IOF_5+IGM_0*IOGF_5+5*IGM_1*IOF_4+5*IGM_1*IOGF_4+10*IGM_2*IOF_3+10*IGM_2*IOGF_3+10*IGM_3*IOF_2+10*IGM_3*IOGF_2+5*IGM_4*IOF_1+5*IGM_4*IOGF_1)*betaOOFM^2+(gammaGM^2+mu^2+nuGOM^2)*IGM_5-gammaOM^2*IOGM_5+IGM_6, ((-IGF_0-IOGF_0)*IOM_5-IOGF_5*IOM_0-IGF_5*IOM_0-5*IOGF_4*IOM_1-5*IGF_4*IOM_1-10*IOGF_3*IOM_2-10*IGF_3*IOM_2-10*IOGF_2*IOM_3-10*IGF_2*IOM_3-5*IOGF_1*IOM_4-5*IGF_1*IOM_4)*betaGGFM^2+((-IGF_0-IOGF_0)*IGM_5-10*IGM_2*IOGF_3-10*IGF_3*IGM_2-10*IGM_3*IOGF_2-10*IGF_2*IGM_3-5*IGM_4*IOGF_1-5*IGF_1*IGM_4-5*IGM_1*IOGF_4-IGM_0*IOGF_5-5*IGF_4*IGM_1-IGF_5*IGM_0)*betaGOFM^2+((-IOF_0-IOGF_0)*IOM_5-5*IOF_1*IOM_4-10*IOF_2*IOM_3-10*IOF_3*IOM_2-5*IOF_4*IOM_1-IOF_5*IOM_0-IOGF_5*IOM_0-5*IOGF_4*IOM_1-10*IOGF_3*IOM_2-10*IOGF_2*IOM_3-5*IOGF_1*IOM_4)*betaOGFM^2+((-IOF_0-IOGF_0)*IGM_5-10*IGM_2*IOF_3-10*IGM_2*IOGF_3-10*IGM_3*IOF_2-10*IGM_3*IOGF_2-5*IGM_4*IOF_1-5*IGM_4*IOGF_1-5*IGM_1*IOF_4-IGM_0*IOF_5-5*IGM_1*IOGF_4-IGM_0*IOGF_5)*betaOOFM^2-IGM_5*nuGOM^2-IOM_5*nuOGM^2+(gammaGM^2+gammaOM^2+mu^2)*IOGM_5+IOGM_6, ((2*IOGF_0+2*IGF_0)*SM_4^2+(2*IGF_4+2*IOGF_4)*SM_0^2+(8*IGF_3+8*IOGF_3)*SM_1^2+(12*IGF_2+12*IOGF_2)*SM_2^2+8*SM_3^2*(IOGF_1+IGF_1))*betaGGFM^2+((2*IOGF_0+2*IGF_0)*SM_4^2+(2*IGF_4+2*IOGF_4)*SM_0^2+(8*IGF_3+8*IOGF_3)*SM_1^2+(12*IGF_2+12*IOGF_2)*SM_2^2+8*SM_3^2*(IOGF_1+IGF_1))*betaGOFM^2+((2*IOGF_0+2*IOF_0)*SM_4^2+(2*IOF_4+2*IOGF_4)*SM_0^2+(8*IOF_3+8*IOGF_3)*SM_1^2+(12*IOF_2+12*IOGF_2)*SM_2^2+8*SM_3^2*(IOGF_1+IOF_1))*betaOGFM^2+((2*IOGF_0+2*IOF_0)*SM_4^2+(2*IOF_4+2*IOGF_4)*SM_0^2+(8*IOF_3+8*IOGF_3)*SM_1^2+(12*IOF_2+12*IOGF_2)*SM_2^2+8*SM_3^2*(IOGF_1+IOF_1))*betaOOFM^2+2*SM_4^2*mu^2-2*gammaGM^2*IGM_4-2*gammaOM^2*IOM_4+2*SM_5^2, -IOGM_5-IOM_5+115942606718361565920408055267237805352214355162586072833382184303857490913393898218777944130977113148035919454019044410945470607405121829282274007707409514029754698145935429107786458517749441200039305162529836207794666164012335665605237455188016018755171128653948812451225860646961451367412282446678277542706188180961474083/4, ((IGF_0+IOGF_0)*IOM_5+10*IGF_2*IOM_3+10*IGF_3*IOM_2+5*IGF_4*IOM_1+IGF_5*IOM_0+IOGF_5*IOM_0+5*IOGF_4*IOM_1+10*IOGF_3*IOM_2+10*IOGF_2*IOM_3+5*IOGF_1*IOM_4+5*IGF_1*IOM_4)*betaGGFM^2+(-IGF_0*SM_5^2-5*IGF_1*SM_4^2-10*IGF_2*SM_3^2-10*IGF_3*SM_2^2-5*IGF_4*SM_1^2-IGF_5*SM_0^2-IOGF_0*SM_5^2-5*IOGF_1*SM_4^2-10*IOGF_2*SM_3^2-10*IOGF_3*SM_2^2-5*IOGF_4*SM_1^2-IOGF_5*SM_0^2)*betaGOFM^2+((IOF_0+IOGF_0)*IOM_5+5*IOF_1*IOM_4+10*IOF_2*IOM_3+10*IOF_3*IOM_2+5*IOF_4*IOM_1+IOF_5*IOM_0+IOGF_5*IOM_0+5*IOGF_4*IOM_1+10*IOGF_3*IOM_2+10*IOGF_2*IOM_3+5*IOGF_1*IOM_4)*betaOGFM^2+(-IOF_0*SM_5^2-5*IOF_1*SM_4^2-10*IOF_2*SM_3^2-10*IOF_3*SM_2^2-5*IOF_4*SM_1^2-IOF_5*SM_0^2-IOGF_0*SM_5^2-5*IOGF_1*SM_4^2-10*IOGF_2*SM_3^2-10*IOGF_3*SM_2^2-5*IOGF_4*SM_1^2-IOGF_5*SM_0^2)*betaOOFM^2+(gammaOM^2+mu^2+nuOGM^2)*IOM_5-gammaGM^2*IOGM_5+IOM_6, -IGF_5-IOGF_5+243194413437837430527356025409270315671447903004069529579152707776652275783719423671687889942854827725497214914185302006066241550176759557835026882836186000171094329710348914795255156051109870664389128041370035859785206260465315043174331188653921323205197688203195719007034607154418626115600250684542559200338880631930732037/4, (-IGM_0*SF_5^2-5*IGM_1*SF_4^2-10*IGM_2*SF_3^2-10*IGM_3*SF_2^2-5*IGM_4*SF_1^2-IGM_5*SF_0^2-IOGM_0*SF_5^2-5*IOGM_1*SF_4^2-10*IOGM_2*SF_3^2-10*IOGM_3*SF_2^2-5*IOGM_4*SF_1^2-IOGM_5*SF_0^2)*betaGGMF^2+((IGM_0+IOGM_0)*IGF_5+10*IGF_2*IGM_3+10*IGF_2*IOGM_3+10*IGF_3*IGM_2+10*IGF_3*IOGM_2+5*IGF_4*IGM_1+5*IGF_4*IOGM_1+5*IGF_1*IGM_4+IGF_0*IGM_5+5*IGF_1*IOGM_4+IGF_0*IOGM_5)*betaGOMF^2+(-IOGM_0*SF_5^2-5*IOGM_1*SF_4^2-10*IOGM_2*SF_3^2-10*IOGM_3*SF_2^2-5*IOGM_4*SF_1^2-IOGM_5*SF_0^2-IOM_0*SF_5^2-5*IOM_1*SF_4^2-10*IOM_2*SF_3^2-10*IOM_3*SF_2^2-5*IOM_4*SF_1^2-IOM_5*SF_0^2)*betaOGMF^2+((IOM_0+IOGM_0)*IGF_5+10*IGF_2*IOM_3+10*IGF_2*IOGM_3+10*IGF_3*IOM_2+10*IGF_3*IOGM_2+5*IGF_4*IOM_1+5*IGF_4*IOGM_1+5*IGF_1*IOM_4+IOM_5*IGF_0+5*IGF_1*IOGM_4+IGF_0*IOGM_5)*betaOOMF^2+(gammaGF^2+mu^2+nuGOF^2)*IGF_5-gammaOF^2*IOGF_5+IGF_6, ((-IGM_0-IOGM_0)*IOF_5-10*IGM_2*IOF_3-10*IGM_3*IOF_2-5*IGM_4*IOF_1-IGM_5*IOF_0-IOF_0*IOGM_5-5*IOF_1*IOGM_4-10*IOF_2*IOGM_3-10*IOF_3*IOGM_2-5*IOF_4*IOGM_1-5*IGM_1*IOF_4)*betaGGMF^2+((-IGM_0-IOGM_0)*IGF_5-10*IGF_3*IGM_2-10*IGF_2*IGM_3-5*IGF_1*IGM_4-IGF_0*IGM_5-IGF_0*IOGM_5-5*IGF_1*IOGM_4-5*IGF_4*IOGM_1-10*IGF_3*IOGM_2-10*IGF_2*IOGM_3-5*IGF_4*IGM_1)*betaGOMF^2+((-IOGM_0-IOM_0)*IOF_5-IOF_0*IOM_5-IOF_0*IOGM_5-5*IOF_1*IOM_4-5*IOF_1*IOGM_4-10*IOF_2*IOM_3-10*IOF_2*IOGM_3-10*IOF_3*IOM_2-10*IOF_3*IOGM_2-5*IOF_4*IOM_1-5*IOF_4*IOGM_1)*betaOGMF^2+((-IOGM_0-IOM_0)*IGF_5-IOM_5*IGF_0-IGF_0*IOGM_5-5*IGF_4*IOM_1-10*IGF_3*IOM_2-10*IGF_2*IOM_3-5*IGF_1*IOM_4-5*IGF_1*IOGM_4-5*IGF_4*IOGM_1-10*IGF_3*IOGM_2-10*IGF_2*IOGM_3)*betaOOMF^2-IGF_5*nuGOF^2-IOF_5*nuOGF^2+(gammaGF^2+gammaOF^2+mu^2)*IOGF_5+IOGF_6, ((2*IOGM_0+2*IGM_0)*SF_4^2+(2*IGM_4+2*IOGM_4)*SF_0^2+(8*IGM_3+8*IOGM_3)*SF_1^2+(12*IGM_2+12*IOGM_2)*SF_2^2+8*SF_3^2*(IOGM_1+IGM_1))*betaGGMF^2+((2*IOGM_0+2*IGM_0)*SF_4^2+(2*IGM_4+2*IOGM_4)*SF_0^2+(8*IGM_3+8*IOGM_3)*SF_1^2+(12*IGM_2+12*IOGM_2)*SF_2^2+8*SF_3^2*(IOGM_1+IGM_1))*betaGOMF^2+((2*IOGM_0+2*IOM_0)*SF_4^2+(2*IOM_4+2*IOGM_4)*SF_0^2+(8*IOM_3+8*IOGM_3)*SF_1^2+(12*IOM_2+12*IOGM_2)*SF_2^2+8*SF_3^2*(IOM_1+IOGM_1))*betaOGMF^2+((2*IOGM_0+2*IOM_0)*SF_4^2+(2*IOM_4+2*IOGM_4)*SF_0^2+(8*IOM_3+8*IOGM_3)*SF_1^2+(12*IOM_2+12*IOGM_2)*SF_2^2+8*SF_3^2*(IOM_1+IOGM_1))*betaOOMF^2+2*SF_4^2*mu^2-2*gammaOF^2*IOF_4-2*gammaGF^2*IGF_4+2*SF_5^2, -IOF_5-IOGF_5+53146577431679143305609103693626059218886129811989031767471873462410155486313300494896011231400205146770932302308695049095411947696075932175397710967967879056608464891997345598672287148268824148722202715101765055981636213329819480903398351739822878400331302449987855652981976796491468290239705790822900682816070441922471713/4, ((IGM_0+IOGM_0)*IOF_5+10*IGM_2*IOF_3+10*IGM_3*IOF_2+5*IGM_4*IOF_1+IGM_5*IOF_0+IOF_0*IOGM_5+5*IOF_1*IOGM_4+10*IOF_2*IOGM_3+10*IOF_3*IOGM_2+5*IOF_4*IOGM_1+5*IGM_1*IOF_4)*betaGGMF^2+(-IGM_0*SF_5^2-5*IGM_1*SF_4^2-10*IGM_2*SF_3^2-10*IGM_3*SF_2^2-5*IGM_4*SF_1^2-IGM_5*SF_0^2-IOGM_0*SF_5^2-5*IOGM_1*SF_4^2-10*IOGM_2*SF_3^2-10*IOGM_3*SF_2^2-5*IOGM_4*SF_1^2-IOGM_5*SF_0^2)*betaGOMF^2+((IOM_0+IOGM_0)*IOF_5+IOF_0*IOM_5+IOF_0*IOGM_5+5*IOF_1*IOM_4+5*IOF_1*IOGM_4+10*IOF_2*IOM_3+10*IOF_2*IOGM_3+10*IOF_3*IOM_2+10*IOF_3*IOGM_2+5*IOF_4*IOM_1+5*IOF_4*IOGM_1)*betaOGMF^2+(-IOGM_0*SF_5^2-5*IOGM_1*SF_4^2-10*IOGM_2*SF_3^2-10*IOGM_3*SF_2^2-5*IOGM_4*SF_1^2-IOGM_5*SF_0^2-IOM_0*SF_5^2-5*IOM_1*SF_4^2-10*IOM_2*SF_3^2-10*IOM_3*SF_2^2-5*IOM_4*SF_1^2-IOM_5*SF_0^2)*betaOOMF^2+(gammaOF^2+mu^2+nuOGF^2)*IOF_5-gammaGF^2*IOGF_5+IOF_6, -IGM_6-IOGM_6-73987144159660152478926141047843834454812749398180698314580214996910313313332136628267813125055157904138275740635412275237303633150279600825529616450298407856896428649055972865227899234401468787937956104117174976535725584491375909690392012702991969607170828929906540518861745344071918819419452524451769988102725138989858210527406993850483385419871150193559940374957503396634670354209/4, (-IGF_0*SM_6^2-6*IGF_1*SM_5^2-15*IGF_2*SM_4^2-20*IGF_3*SM_3^2-15*IGF_4*SM_2^2-6*IGF_5*SM_1^2-IGF_6*SM_0^2-IOGF_0*SM_6^2-6*IOGF_1*SM_5^2-15*IOGF_2*SM_4^2-20*IOGF_3*SM_3^2-15*IOGF_4*SM_2^2-6*IOGF_5*SM_1^2-IOGF_6*SM_0^2)*betaGGFM^2+((IGF_0+IOGF_0)*IGM_6+6*IGF_1*IGM_5+15*IGF_2*IGM_4+20*IGF_3*IGM_3+15*IGF_4*IGM_2+6*IGF_5*IGM_1+IGF_6*IGM_0+IGM_0*IOGF_6+6*IGM_1*IOGF_5+15*IGM_2*IOGF_4+20*IOGF_3*IGM_3+15*IGM_4*IOGF_2+6*IGM_5*IOGF_1)*betaGOFM^2+(-IOF_0*SM_6^2-6*IOF_1*SM_5^2-15*IOF_2*SM_4^2-20*IOF_3*SM_3^2-15*IOF_4*SM_2^2-6*IOF_5*SM_1^2-IOF_6*SM_0^2-IOGF_0*SM_6^2-6*IOGF_1*SM_5^2-15*IOGF_2*SM_4^2-20*IOGF_3*SM_3^2-15*IOGF_4*SM_2^2-6*IOGF_5*SM_1^2-IOGF_6*SM_0^2)*betaOGFM^2+((IOF_0+IOGF_0)*IGM_6+IOF_6*IGM_0+IGM_0*IOGF_6+6*IGM_1*IOF_5+6*IGM_1*IOGF_5+15*IGM_2*IOF_4+15*IGM_2*IOGF_4+20*IGM_3*IOF_3+20*IOGF_3*IGM_3+15*IGM_4*IOF_2+15*IGM_4*IOGF_2+6*IGM_5*IOF_1+6*IGM_5*IOGF_1)*betaOOFM^2+(gammaGM^2+mu^2+nuGOM^2)*IGM_6-gammaOM^2*IOGM_6+IGM_7, ((-IGF_0-IOGF_0)*IOM_6-IOM_0*IOGF_6-IGF_6*IOM_0-6*IOM_1*IOGF_5-6*IGF_5*IOM_1-15*IOGF_4*IOM_2-15*IGF_4*IOM_2-20*IOGF_3*IOM_3-20*IGF_3*IOM_3-15*IOGF_2*IOM_4-15*IOM_4*IGF_2-6*IOM_5*IGF_1-6*IOGF_1*IOM_5)*betaGGFM^2+((-IGF_0-IOGF_0)*IGM_6-6*IGM_1*IOGF_5-6*IGF_5*IGM_1-15*IGM_2*IOGF_4-15*IGF_4*IGM_2-20*IOGF_3*IGM_3-20*IGF_3*IGM_3-15*IGM_4*IOGF_2-15*IGF_2*IGM_4-6*IGF_1*IGM_5-6*IGM_5*IOGF_1-IGM_0*IOGF_6-IGF_6*IGM_0)*betaGOFM^2+((-IOF_0-IOGF_0)*IOM_6-6*IOM_5*IOF_1-15*IOM_4*IOF_2-20*IOM_3*IOF_3-15*IOM_2*IOF_4-6*IOM_1*IOF_5-IOF_6*IOM_0-IOM_0*IOGF_6-6*IOM_1*IOGF_5-15*IOGF_4*IOM_2-20*IOGF_3*IOM_3-15*IOGF_2*IOM_4-6*IOGF_1*IOM_5)*betaOGFM^2+((-IOF_0-IOGF_0)*IGM_6-6*IGM_1*IOF_5-6*IGM_1*IOGF_5-15*IGM_2*IOF_4-15*IGM_2*IOGF_4-20*IGM_3*IOF_3-20*IOGF_3*IGM_3-15*IGM_4*IOF_2-15*IGM_4*IOGF_2-6*IGM_5*IOF_1-6*IGM_5*IOGF_1-IOF_6*IGM_0-IGM_0*IOGF_6)*betaOOFM^2-IGM_6*nuGOM^2-IOM_6*nuOGM^2+(gammaGM^2+gammaOM^2+mu^2)*IOGM_6+IOGM_7, ((2*IOGF_0+2*IGF_0)*SM_5^2+(2*IGF_5+2*IOGF_5)*SM_0^2+(10*IGF_4+10*IOGF_4)*SM_1^2+(20*IGF_3+20*IOGF_3)*SM_2^2+(20*IGF_2+20*IOGF_2)*SM_3^2+10*SM_4^2*(IOGF_1+IGF_1))*betaGGFM^2+((2*IOGF_0+2*IGF_0)*SM_5^2+(2*IGF_5+2*IOGF_5)*SM_0^2+(10*IGF_4+10*IOGF_4)*SM_1^2+(20*IGF_3+20*IOGF_3)*SM_2^2+(20*IGF_2+20*IOGF_2)*SM_3^2+10*SM_4^2*(IOGF_1+IGF_1))*betaGOFM^2+((2*IOGF_0+2*IOF_0)*SM_5^2+(2*IOF_5+2*IOGF_5)*SM_0^2+(10*IOF_4+10*IOGF_4)*SM_1^2+(20*IOF_3+20*IOGF_3)*SM_2^2+(20*IOF_2+20*IOGF_2)*SM_3^2+10*SM_4^2*(IOGF_1+IOF_1))*betaOGFM^2+((2*IOGF_0+2*IOF_0)*SM_5^2+(2*IOF_5+2*IOGF_5)*SM_0^2+(10*IOF_4+10*IOGF_4)*SM_1^2+(20*IOF_3+20*IOGF_3)*SM_2^2+(20*IOF_2+20*IOGF_2)*SM_3^2+10*SM_4^2*(IOGF_1+IOF_1))*betaOOFM^2+2*SM_5^2*mu^2-2*gammaGM^2*IGM_5-2*gammaOM^2*IOM_5+2*SM_6^2, -IGM_7-IOGM_7+12278367466658149075069043646393160189618471804170622262029361875696109157111506078563836604090766770141501490593553424628966767576991118091594423292920921187853034215252960619747868068019251394215497452037734733851618430867633557871777210441159202691136340938672956206303024904170716408836823832777449189211068645657460268125265339986817454510327690824477973690271262676072417823156209963616712194113018311783855733186835743429548810542541853/8, -IOGM_6-IOM_6-26196543006826920580137742534939701216354646361048874228186620341887680444468213754479701266387371678523429792046479076633597329616125240412807631725789286917306104218623647783999783511530468406946240383089560746549287810715736291311359894759255834352603199452152367970560629457070651771952377726177447038222732824382132793917807095152787757818678224604237884909472355816449562462353/8, -IGF_6-IOGF_6-94584984726043836124121819725752954416396611323406198487090413226091851237469150509119089441953459598610027834578510038233846055568066091468863009481190842150920426629382394699617676588262564802478749655027933540392353653861058006897661614118109490795065689264044297022749234586781342949157356445563893760520574447659513839379377350999958376250279743281929006730491583287227117375687/8, -IOF_6-IOGF_6-10543581345582237533534572730041841653388210125743530749985197403994051907723955710797889833388342434309814238000422016654466016341128763299955983012826500836725892211072445126169498652190077663755457160438897905516267385628182688034685940179000662747651645180926802074952382076209964831733473152100996939134833521703958726602378971356832723590241484802778515361698147010124834780989/4, 2*z_aux-1];
vars:=[IOGM_7, IGM_7, SM_6, IOM_6, IOGM_6, IOGF_6, IOF_6, IGM_6, IGF_6, SM_5, SF_5, IOM_5, IOGM_5, IOGF_5, IOF_5, IGM_5, IGF_5, SM_4, SF_4, IOM_4, IOGM_4, IOGF_4, IOF_4, IGM_4, IGF_4, SM_3, SF_3, IOM_3, IOGM_3, IOGF_3, IOF_3, IGM_3, IGF_3, SM_2, SF_2, IOM_2, IOGM_2, IOGF_2, IOF_2, IGM_2, IGF_2, SM_1, SF_1, IOM_1, IOGM_1, IOGF_1, IOF_1, IGM_1, IGF_1, SM_0, SF_0, IOM_0, IOGM_0, IOGF_0, IOF_0, IGM_0, IGF_0, z_aux, w_aux, gammaGF, gammaGM, gammaOF, gammaOM, mu, nuGOF, nuGOM, nuOGF, nuOGM, betaGGFM, betaGGMF, betaGOFM, betaGOMF, betaOGFM, betaOGMF, betaOOFM, betaOOMF];
gb:=CodeTools[Usage](Groebner[Basis](et_hat, tdeg(op(vars)), characteristic=11863279),output='all');
# [z_aux = z_aux, IOF_6 = IOF_6, IOM_3 = IOM_3, IGF_0 = IGF_0, IGM_4 = IGM_4, IOGF_5 = IOGF_5, gammaOF = gammaOF^2, SF_4 = SF_4^2, betaOOMF = betaOOMF^2, IOGF_1 = IOGF_1, SM_5 = SM_5^2, IGF_3 = IGF_3, IGM_3 = IGM_3, IOGF_4 = IOGF_4, SF_3 = SF_3^2, nuOGF = nuOGF^2, IGF_1 = IGF_1, betaOGFM = betaOGFM^2, SM_0 = SM_0^2, mu = mu^2, SF_5 = SF_5^2, IOM_5 = IOM_5, IOF_1 = IOF_1, IGM_0 = IGM_0, nuGOF = nuGOF^2, betaGOMF = betaGOMF^2, IOF_5 = IOF_5, SM_2 = SM_2^2, SF_1 = SF_1^2, IGF_5 = IGF_5, IOGF_3 = IOGF_3, IOF_2 = IOF_2, IOGF_0 = IOGF_0, IOGM_4 = IOGM_4, IOGM_3 = IOGM_3, betaGGMF = betaGGMF^2, betaOOFM = betaOOFM^2, gammaOM = gammaOM^2, IGF_2 = IGF_2, IOGM_1 = IOGM_1, gammaGF = gammaGF^2, IOGM_5 = IOGM_5, IGM_6 = IGM_6, IGF_4 = IGF_4, IGM_5 = IGM_5, IOGF_2 = IOGF_2, IGM_1 = IGM_1, IOM_4 = IOM_4, SF_0 = SF_0^2, SM_4 = SM_4^2, SM_3 = SM_3^2, betaGGFM = betaGGFM^2, IOGM_6 = IOGM_6, IOM_2 = IOM_2, IOF_4 = IOF_4, betaOGMF = betaOGMF^2, gammaGM = gammaGM^2, IGF_6 = IGF_6, IOM_1 = IOM_1, nuGOM = nuGOM^2, IOM_0 = IOM_0, IOGM_2 = IOGM_2, IOGF_6 = IOGF_6, IOGM_7 = IOGM_7, SM_1 = SM_1^2, nuOGM = nuOGM^2, IOGM_0 = IOGM_0, IOM_6 = IOM_6, SM_6 = SM_6^2, IGM_2 = IGM_2, IOF_3 = IOF_3, betaGOFM = betaGOFM^2, SF_2 = SF_2^2, IGM_7 = IGM_7, IOF_0 = IOF_0]
quit;