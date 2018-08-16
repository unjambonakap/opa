function main ()
  n=100;
  bias=0.9^60;
  p_success=bias/2+0.5;
  p=1-p_success;

  var_x=p*(1-p);

  u=(n/2-p*n)/sqrt(n*var_x);

  prob=1/2*erfc(u);
  prob


  %test2();
  test3();
end

function test2()
  bias=0.4^4
  p1=bias/2+0.5;
  p1
  p2=0.5;
  nc=1000;

  u1=p1*nc;
  v1=nc*p1*(1-p1);

  u2=p2*nc;
  v2=nc*p2*(1-p2);
  u2=(u2+u1)/2;
  v2=(v2+v1)/2;


  u12=u1-u2;
  v12=v1+v2;
  u12
  v12
  pfail=1/2*(erfc(-(0-u12)/sqrt(2*v12)));
  n=2^27;

  pfail
  res=binocdf(1, n, pfail);
  res
end

function [u0,v0,u1,v1]=params(n, bias)
  p0=(1+bias)/2;
  p0
  p1=1/2;

  u0=n*p0;
  u1=n*p1;
  v0=p0*(1-p0)*n;
  v1=p1*(1-p1)*n;

end

function test3()
  bias=0.6^4;
  n=200*8*20;
  n=hex2dec('1102')
  bias=0.64*2-1
  [u0,v0,u1,v1]=params(n, bias);

  m=2^28;

  c=[1/2/v1-1/2/v0, u0/v0-u1/v1, u1*u1/2/v1-u0*u0/2/v0+1/2*log(v1/v0)-log(m-1)];
  r=roots(c);
  r
  u0
  u1
  res=min(r);
  res
  pc=1/2*erfc((res-u0)/sqrt(v0*2));
  pc
  pc2=1/2*erfc((res-u1)/sqrt(v1*2));
  pc2

  %tot=hex2dec('1fd');
  %obs0=hex2dec('173');
  %obs1=hex2dec('119');



  %n0=(obs0+tot)/2;
  %n1=(obs1+tot)/2;

  %obs0/tot
  %obs1/tot
  %bias


  %[u0,v0,u1,v1]=params(tot, bias);
  %u0
  %u1
  %n0
  %n1
  %v0
  %v1
  %tot

  %prob1=1/2*erfc((n1-u1)/sqrt(2*v1));
  %prob0=1/2*erfc((n0-u0)/sqrt(2*v0));
  %prob1
  %prob0


end
