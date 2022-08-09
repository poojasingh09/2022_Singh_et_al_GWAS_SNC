##pooja.singh09@gmail.com
##may 2021
##figure2A venn diagram was made in Venny https://bioinfogp.cnb.csic.es/tools/venny/ and expected values and overlap significance was calcualted using phyper code below




#phyper gene overlap


#SNC TR IR

phyper(76, 436+76+2, 31954, 534+76+4, lower.tail = F, log.p = FALSE)
[1] 1.025901e-45

((436+76+2)*(534+76+4))/31954
[1] 9.876573

#SNC IT IR

phyper(4, 76+4+2, 31954, 534+76+4, lower.tail = F, log.p = FALSE)
[1] 0.02076045

((76+4+2)*(534+76+4))/31954
[1] 1.57564


#SNC IT TR

phyper(2, 76+4+2, 31954, 436+76+2, lower.tail = F, log.p = FALSE)
[1] 0.1449994

((76+4+2)*(436+76+2))/31954
[1] 1.319021


