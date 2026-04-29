function sudan(x,y,n) {
  if (n == 1) {
    return OmegaNum.add(OmegaNum.mul(OmegaNum.pow(2,y), OmegaNum.add(x,2)),OmegaNum.sub(y,2));
  } else if (OmegaNum.cmp(y,0) == 1) {
    if (OmegaNum.cmp(6,y) == 1) {
      let z = sudan(x,OmegaNum.sub(y,1),n);
      return sudan(z,OmegaNum.add(z,OmegaNum.add(y,1)),n-1);
    } else {
      if (n == 2) {
        return OmegaNum.tetr(2,OmegaNum.sub(y,4),sudan(x,4,2));
      } else {
        return OmegaNum.arrow(sudan(x,2,n),n,OmegaNum.sub(y,1));
      }
    } 
  } else {
    return x;
  }
}

return sudan(7,3,2);  //replace this by any other value