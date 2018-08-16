#ifndef _H_GFQ_MAPPER
#define _H_GFQ_MAPPER

#define GFQ_TEST 1

//void GF_qZech::construct(const Field<int> &baseField, const Poly &p){
//
//    int qm=getSize();
//    PolyRing pr(baseField);
//
//    Poly w=utils::findPrimitiveElem(baseField, p);
//
//    m_qm1=qm-1;
//    zechslog=(int*)malloc((m_qm1+1)*sizeof(int));
//    lst_int=(int*)malloc(m_car*sizeof(int));
//
//
//    std::map<Poly,int> mp;
//    std::map<int,Poly> test;
//#if GFQ_TEST
//    test[0]=Poly::zero();
//    mp[Poly::zero()]=0;
//#endif
//
//    Poly cur=Poly::one();
//    for (int i=0; i<qm-1; ++i){
//        mp[cur]=i+1;
//        printf("%d >>> ", i+1); cur.disp();
//#if GFQ_TEST
//        test[i+1]=cur;
//#endif
//        cur=pr.mulmod(cur,w, p);
//    }
//
//    std::map<Poly,int>::iterator it;
//    for (it=mp.begin(); it!=mp.end(); ++it){
//        Poly u=pr.add(it->first, Poly::one());
//
//        assert(mp.count(u));
//        zechslog[it->second]=mp[u];
//    }
//
//    for (int i=0; i<m_car; ++i){
//        lst_int[i]=mp[pr.getu32(i)];
//    }
//
//#if GFQ_TEST
//    for (int i=1; i<qm; ++i){
//        for (int j=1; j<qm; ++j){
//            Poly res=pr.mulmod(test[i], test[j], p);
//            //printf("======\n");
//            //test[i].disp();
//            //test[j].disp();
//            //res.disp();
//            int ans=mul(i,j);
//            assert(mp.count(res));
//            int sol=mp[res];
//
//            //printf(" %d\n", zechslog[0]);
//            //printf("GET >> %d %d\n",ans,sol);
//            assert(ans==sol);
//
//        }
//    }
//#endif
//
//}


#endif
