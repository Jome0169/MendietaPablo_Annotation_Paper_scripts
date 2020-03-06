 find . -iname "*root*leaf*" -exec rm {} \;
 find . -iname "*leaf*root*" -exec rm {} \;
 find . -iname "*leaf*ear*" -exec rm {} \;
 find . -iname "ear*leaf*" -exec rm {} \;
 find . -iname "*ear*root*" -exec rm {} \;
 find . -iname "*root*ear*" -exec rm {} \;
