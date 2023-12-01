
library(Rcpp)
library(dplyr)
library(igraph)
library(Matrix)
sourceCpp("/gdata01/user/xuhe/family_relatedness/simulation-2023-03-22/code/subfunc-2023-03-22XH.cpp")

###### main function #####

SPA_G = function(obj.precond,   # value of function precond()
                 GenoMat,       # column names of SNP ID and row names of Subject IDs are required
                 Cutoff = 2,
                 impute.method = "fixed",
                 missing.cutoff = 0.15,
                 min.maf = 0.0001,
                 G.model = "Add",
                 zeta = 0, # set zeta = - 0.02 for beta and - 0.0002 for tau is better.
                 tol = 1e-5) # set tol = 1e-4 for beta and 1e-5 for tau.
{
  SubjID = as.character(obj.precond$SubjID)
  GenoMat = GenoMat[SubjID, ]
  
  ## check input
  par.list = list(pwd = getwd(),
                  sessionInfo = sessionInfo(),
                  impute.method = impute.method,
                  missing.cutoff = missing.cutoff,
                  min.maf = min.maf,
                  G.model = G.model)
  
  # check_input(obj.null, GenoMat, par.list)
  if(length(obj.precond$TwoSubj_list) == 0 & length(obj.precond$ThreeSubj_list) == 0)
    warning("There is no relatedness in the data!")
  
  print(paste0("Sample size is ", nrow(GenoMat), "."))
  print(paste0("Number of variants is ", ncol(GenoMat), "."))
  
  ### Prepare the main output data frame
  n.Geno = ncol(GenoMat)
  output = matrix(NA, nrow = n.Geno, ncol = 5)
  colnames(output) = c("MAF", "missing.rate","Score", 
                       "p.value.spa.G.GRM", "p.value.norm.GRM")
  SNPID = colnames(GenoMat)
  
  ### Start analysis
  print("Start Analyzing...")
  print(Sys.time())
  
  # Cycle for genotype matrix
  for(i in 1:n.Geno){
    geno = GenoMat[,i]
    
    output.one.SNP = SPA_G.one.SNP(obj.precond = obj.precond,
                                   geno = geno,
                                   Cutoff = Cutoff,
                                   impute.method = impute.method,
                                   missing.cutoff = missing.cutoff,
                                   min.maf = min.maf,
                                   G.model = G.model,
                                   zeta = zeta,
                                   tol = tol)
    output[i,] = output.one.SNP
    
    if(i %% 1000 == 0) cat("Complete analyzing", i, "SNPs.\n")
  }
  output = cbind(SNPID, output)
  
  print("Analysis Complete.")
  print(Sys.time())
  
  return(output)
}

###### preconditioning part #####
precond = function(GRM.file,       # three columns: "ID1", "ID2", and "Value"
                   GenoMat,        # column names of SNP ID and row names of Subject IDs are required
                   ResidMat,       # two columns: "SubjID", and "Resid". The subjects in "ResidMat" should be also in "GRM.file" and "GenoMat"
                   MAFvec,         # a MAF vector corresponding to GenoMat, so that only related samples' GenoMat is inputted.
                   MaxQuantile = 0.75,
                   MinQuantile = 0.25,
                   OutlierRatio = 1.5,
                   MaxNuminFam = 5,
                   MAF_interval = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5))
{
  print(Sys.time())
  
  GRM = data.table::fread(GRM.file)
  GRM$ID1 = as.character(GRM$ID1); GRM$ID2 = as.character(GRM$ID2)
  
  SubjID.In.Resid = ResidMat$SubjID
  SubjID.In.GRM = unique(c(GRM$ID1, GRM$ID2))
  
  if(any(!SubjID.In.Resid %in% SubjID.In.GRM))
    stop("At least one subject in residual matrix does not have GRM information.")
  
  SubjID = SubjID.In.Resid
  GRM = GRM %>% filter(ID1 %in% SubjID & ID2 %in% SubjID)
  
  # Use residual information to define outliers / non-outliers
  Resid = ResidMat$Resid
  Quant = quantile(Resid, probs = c(MinQuantile, MaxQuantile))
  Range = max(Quant) - min(Quant)
  cutoffVec = c(min(Quant) - OutlierRatio * Range, max(Quant) + OutlierRatio * Range)
  
  cat("cutoffVec:\t",cutoffVec,"\n")
  ResidMat$Outlier = ifelse(Resid < cutoffVec[1] | Resid > cutoffVec[2],
                            TRUE, FALSE)
  
  cat("Outliers information is as below\n")
  print(ResidMat %>% filter(Outlier == TRUE) %>% dplyr::select(SubjID, Resid, Outlier) %>% arrange(Resid))
  
  # Decompose the subjects based on family structure and use a greedy algorithm to reduce family size if needed
  GRM1 = GRM
  GRM1$pos1 = ResidMat$Resid[match(GRM$ID1, ResidMat$SubjID)]
  GRM1$pos2 = ResidMat$Resid[match(GRM$ID2, ResidMat$SubjID)]
  GRM1 = GRM1 %>% mutate(Cov = abs(Value * pos1 * pos2))
  
  edges = t(GRM1[, c("ID1", "ID2")])
  graph_GRM = make_graph(edges, directed = F)
  graph_list_all = graph_GRM %>% decompose()
  graph_length = lapply(graph_list_all, length)
  
  graph_list_1 = graph_list_all[graph_length == 1]
  SubjID.unrelated = lapply(graph_list_1, get.vertex.attribute) %>% unlist(use.names = FALSE)
  ResidMat.unrelated = ResidMat %>% filter(SubjID %in% SubjID.unrelated)
  SubjID.unrelated.nonOutlier = ResidMat.unrelated %>% filter(Outlier == FALSE) %>% select(SubjID) %>% unlist(use.names = F)
  
  # Values used in association analysys
  R_GRM_R = GRM1 %>% filter(ID1 %in% SubjID.unrelated) %>% select(Cov) %>% sum
  sum_R_nonOutlier = ResidMat.unrelated %>% filter(Outlier == FALSE) %>% select(Resid) %>% sum
  R_GRM_R_nonOutlier = GRM1 %>% filter(ID1 %in% SubjID.unrelated.nonOutlier) %>% select(Cov) %>% sum
  Resid.unrelated.outliers = ResidMat.unrelated %>% filter(Outlier == TRUE) %>% select(Resid) %>% unlist(use.names = F)
  R_GRM_R_TwoSubjOutlier = 0; TwoSubj_list = ThreeSubj_list = list();
  
  graph_list_updated = list()
  graph_list = graph_list_all[graph_length > 1]
  nGraph = length(graph_list)
  index.outlier = 1
  
  if(nGraph != 0)
  {
    cat("Start process the related residual information.\n")
    
    for(i in 1:nGraph)
    {
      if(i %% 1000 == 0)
        cat("Processing the related residual information:\t", i,"/",nGraph,"\n")
      
      comp1 = graph_list[[i]]
      comp3 = V(comp1)$name
      
      # Step 0: calculate variance for the family
      pos1 = match(comp3, SubjID.In.Resid)
      outlierInFam = any(ResidMat$Outlier[pos1])
      
      block_GRM = make.block.GRM(comp1, GRM)
      
      R_GRM_R.temp = as.numeric(t(Resid[pos1]) %*% block_GRM %*% Resid[pos1])
      R_GRM_R = R_GRM_R + R_GRM_R.temp
      
      if(!outlierInFam)
      {
        sum_R_nonOutlier = sum_R_nonOutlier + sum(ResidMat$Resid[pos1])
        R_GRM_R_nonOutlier = R_GRM_R_nonOutlier + R_GRM_R.temp
        next
      }
      
      # cat("Family ", i, " (with outliers) includes ", length(comp3), " subjects:", comp3, "\n")
      
      vcount = vcount(comp1)   # number of vertices 
      
      if(vcount <= MaxNuminFam)
      {
        graph_list_updated[[index.outlier]] = comp1
        index.outlier = index.outlier + 1
        next
      }
      
      # Step 1: remove the edges until the largest family size is <= MaxNuminFam, default is 5.
      
      comp1.temp = comp1
      GRM1.temp = GRM1 %>% filter(ID1 %in% comp3 | ID2 %in% comp3) %>% arrange(Cov)
      for(j in 1:nrow(GRM1.temp))
      {
        # cat("j:\t",j,"\n")
        edgesToRemove = paste0(GRM1.temp$ID1[j],"|",GRM1.temp$ID2[j])
        comp1.temp = delete.edges(comp1.temp, edgesToRemove)
        vcount = decompose(comp1.temp) %>% sapply(vcount)  # vertices count for the new graph after edge removal
        # cat("vcount:\t",vcount,"\n")
        if(max(vcount) <= MaxNuminFam)
          break;
      }
      
      # cat("Edge removal complete. Counts of vertices:\t", vcount,"\n")
      
      # Step 2: add the (removed) edges while keeping the largest family size <= MaxNuminFam, default is 5.
      
      GRM1.temp = GRM1.temp[1:j,] %>% arrange(desc(Cov))
      comp1 = comp1.temp
      for(k in 1:nrow(GRM1.temp))
      {
        # cat("k:\t",k,"\n")
        edgesToAdd = c(GRM1.temp$ID1[k], GRM1.temp$ID2[k])
        comp1.temp = add.edges(comp1, edgesToAdd)
        
        vcount = decompose(comp1.temp) %>% sapply(vcount)  # vertices count for the new graph after edge removal
        # cat("vcount:\t",vcount,"\n")
        
        if(max(vcount) <= MaxNuminFam)
          comp1 = comp1.temp
      }
      
      comp1 = decompose(comp1)
      
      # cat("Edge add complete. Counts of vertices:\t", comp1 %>% sapply(vcount),"\n")
      
      for(k in 1:length(comp1))
      {
        comp11 = comp1[[k]]
        comp13 = V(comp11)$name
        
        pos2 = match(comp13, SubjID.In.Resid)
        outlierInFam = any(ResidMat$Outlier[pos2])
        
        block_GRM = make.block.GRM(comp11, GRM)
        
        R_GRM_R.temp = as.numeric(t(Resid[pos2]) %*% block_GRM %*% Resid[pos2])
        
        if(!outlierInFam){
          sum_R_nonOutlier = sum_R_nonOutlier + sum(ResidMat$Resid[pos2])
          R_GRM_R_nonOutlier = R_GRM_R_nonOutlier + R_GRM_R.temp
        }else{
          graph_list_updated[[index.outlier]] = comp11
          index.outlier = index.outlier + 1;
        }
      }
    }
    
    cat("Start process the Chow-Liu tree.\n")
    
    # Make a list of array index.
    arr.index = list()
    for(n in 1:MaxNuminFam)
    {
      temp = c()
      for(i in 1:n)
      {
        indexString = rep("c(1, 1, 1)", n)
        indexString[i] = "0:2"
        indexString = paste0(indexString, collapse = "%o%")
        cmd = paste0("temp = c(temp, list(arr.index", i, "=", indexString, "))")
        eval(parse(text = cmd))
      }
      arr.index[[n]] = temp
    }
    
    # build chou-liu-tree.
    n.outliers = length(graph_list_updated)
    if(n.outliers != 0)
    {
      ## The below values are only used in chou.liu.tree
      # MAFvec = colMeans(GenoMat, na.rm = T)/2
      pro.var = 2 * MAFvec^2 * (1 - MAFvec)^2
      wi = sqrt(pro.var/(1 - pro.var))
      
      TwofamID.index = ThreefamID.index = 0
      for(index.outlier in 1:n.outliers)
      {
        if(index.outlier %% 1000 == 0)
          cat("Processing the CLT for families with outliers:\t", TwofamID.index, ", ", ThreefamID.index, "/", nGraph, "\n")
        
        comp1 = graph_list_updated[[index.outlier]]
        comp3 = V(comp1)$name
        n1 = length(comp3)
        pos1 = match(comp3, SubjID.In.Resid)
        
        Resid.temp = ResidMat$Resid[pos1]
        
        if(n1 == 1)
        {
          Resid.unrelated.outliers = c(Resid.unrelated.outliers, Resid.temp)
          next;
        }
        
        block_GRM = make.block.GRM(comp1, GRM)
        
        CLT = chou.liu.tree(GRM = block_GRM,
                            GenoMat = GenoMat[comp3,],
                            pro.var = pro.var,
                            wi = wi,
                            MAF_interval = MAF_interval)
        
        if(n1 == 2)
        {
          TwofamID.index = TwofamID.index + 1;
          
          R_GRM_R_TwoSubjOutlier.temp = as.numeric(t(Resid.temp) %*% block_GRM %*% Resid.temp)
          R_GRM_R_TwoSubjOutlier = R_GRM_R_TwoSubjOutlier + R_GRM_R_TwoSubjOutlier.temp
          
          Rho.temp = CLT$pa + 0.5*CLT$pb
          midterm = sqrt(Rho.temp^2 - CLT$pa)
          
          TwoSubj_list[[TwofamID.index]] = list(Resid = Resid.temp,
                                                Rho = c(Rho.temp + midterm, Rho.temp - midterm))
          next;
        }
        
        ThreefamID.index = ThreefamID.index + 1;
        
        stand.S.temp = array(rowSums(mapply(function(x, y) x*y, arr.index[[n1]], Resid.temp)), rep(3, n1))
        
        ThreeSubj_list[[ThreefamID.index]] = list(CLT = CLT,
                                                  stand.S = stand.S.temp)
      }
      cat("Completed processing the CLT for families with outliers:\t", TwofamID.index, ", ", ThreefamID.index, "/", nGraph, "\n")
    }
  }
  
  print(Sys.time())
  
  obj.precond = list(Resid = Resid, SubjID = SubjID, Resid.unrelated.outliers = Resid.unrelated.outliers,
                     R_GRM_R = R_GRM_R, R_GRM_R_TwoSubjOutlier = R_GRM_R_TwoSubjOutlier,
                     sum_R_nonOutlier = sum_R_nonOutlier, R_GRM_R_nonOutlier = R_GRM_R_nonOutlier,
                     TwoSubj_list = TwoSubj_list, ThreeSubj_list = ThreeSubj_list, 
                     MAF_interval = MAF_interval)
  
  return(obj.precond)
}

make.block.GRM = function(graph, 
                          GRM)    # three columns: "ID1", "ID2", and "Value"
{
  comp2 = get.data.frame(graph)
  
  # igraph gives an unexpected additional loop, which may change the block GRM
  # the below is to remove the additional loop
  comp2 = comp2[!duplicated(comp2),]
  
  comp3 = V(graph)$name
  
  colnames(GRM) = c("to", "from", "Value")
  
  n1 = nrow(comp2)
  comp2 = merge(comp2, GRM)
  n2 = nrow(comp2)
  
  if(n1 != n2)
    stop("Ask Wenjian Bi (wenjianb@pku.edu.cn) to check why 'n1 != n2'.")
  
  block_GRM = sparseMatrix(i = match(comp2$from, comp3),
                           j = match(comp2$to, comp3),
                           x = comp2$Value,
                           symmetric = T)
  return(block_GRM)
}

chou.liu.tree = function(GRM,         # block GRM denotes family relatedness
                         GenoMat,     # corresponding GenoMat
                         pro.var,     # 2*pi^2*(1-pi)^2, where pi comes from Binom(2, pi)
                         wi,
                         MAF_interval)
{
  n = nrow(GRM)
  
  pro.frame = data.table::data.table()
  for(i in 1:(n - 1))
  {
    for(j in (i + 1):n)
    {
      num = GenoMat[i,] - GenoMat[j,]
      
      pc = 0.5 * weighted.mean(((abs(num - 1) + abs(num + 1)- 2)/pro.var), wi, na.rm = T)
      pc = ifelse(pc > (1-GRM[i,j])^2, (1-GRM[i,j])^2-1e-10, ifelse(pc < 1-2*GRM[i,j], 1-2*GRM[i,j], pc))
      
      pb = 2 - 2*pc - 2*GRM[i,j]
      
      pa = 2*GRM[i,j] + pc - 1
      
      pro.frame = bind_rows(pro.frame,
                            data.frame(ID1 = i, ID2 = j, pa = pa, pb = pb, pc = pc))
    }
  }
  
  if(n == 2)
  {
    return(pro.frame)
  }else{
    # MAF_interval = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
    CLT = list()
    
    for(index.CLT in 1:length(MAF_interval))
    {
      mu = MAF_interval[index.CLT]
      
      # p = c(G0, G1, G2)
      p0 = c((1-mu)^2, 2*mu*(1-mu), mu^2)
      
      # p = c(G00, G10, G20, G01, G11, G21, G02, G12, G22)
      pa.allele2 = c((1-mu)^2, 0, 0, 0, 2*mu*(1-mu), 0, 0, 0, mu^2)
      
      pb.allele1 = c((1-mu)^3, mu*(1-mu)^2, 0, mu*(1-mu)^2, mu*(1-mu), mu^2*(1-mu), 0, mu^2*(1-mu), mu^3)
      
      pc.allele0 = c((1-mu)^4, 2*mu*(1-mu)^3, mu^2*(1-mu)^2, 2*mu*(1-mu)^3, 4*mu^2*(1-mu)^2, 
                     2*mu^3*(1-mu), mu^2*(1-mu)^2, 2*mu^3*(1-mu), mu^4)
      
      for(j in 1:nrow(pro.frame))
      {
        pro = pro.frame$pa[j] * pa.allele2 + pro.frame$pb[j] * pb.allele1 + pro.frame$pc[j] * pc.allele0
        
        entropy = sum(pro * log(pro/pc.allele0), na.rm = T)
        pro.frame$entropy[j] = entropy
      }
      
      Max_span_tree = pro.frame %>% graph_from_data_frame(directed = T) %>% 
        mst(weights = - pro.frame$entropy, algorithm = "prim") %>% get.edgelist() %>% 
        data.table::as.data.table() %>% transform(V1 = as.numeric(V1), V2 = as.numeric(V2))
      colnames(Max_span_tree) = c("ID1", "ID2")
      
      mst.pro.frame = merge(Max_span_tree, pro.frame, all.x = T)
      
      arr.prob = array(1, dim = rep(3, n))
      for(i in 1:n)
        dimnames(arr.prob)[[i]] = paste0("ID",i,":",0:2) 
      
      vec = c(mst.pro.frame$ID1, mst.pro.frame$ID2); vec = vec[duplicated(vec)]
      
      for(k in 1:(n - 1))
      {
        pro = mst.pro.frame$pa[k] * pa.allele2 + mst.pro.frame$pb[k] * pb.allele1 + mst.pro.frame$pc[k] * pc.allele0
        
        matrix.prob = matrix(pro, 3, 3)
        matrix.index = as.numeric(c(mst.pro.frame$ID1[k], mst.pro.frame$ID2[k]))
        for(i in 1:3){
          for(j in 1:3){
            indexString = rep("", n)
            indexString[matrix.index[1]] = i
            indexString[matrix.index[2]] = j
            indexString = paste0(indexString, collapse = ",")
            cmd = paste0("arr.prob[",indexString,"] = arr.prob[", indexString, "] * matrix.prob[",i,",",j,"]")
            # "arr.prob[1,1,] = arr.prob[1,1,] * matrix.prob[1,1]"
            eval(parse(text = cmd))
          }
        }
      }
      
      for(k in 1:(n - 2))
      {
        vector.prob = p0
        vector.index = vec[k]
        for(i in 1:3){
          indexString = rep("", n)
          indexString[vector.index] = i
          indexString = paste0(indexString, collapse = ",")
          cmd = paste0("arr.prob[",indexString,"] = arr.prob[", indexString, "] / vector.prob[",i,"]")
          # arr.prob[,,1] = arr.prob[,,1] / vector.prob[1]"
          eval(parse(text = cmd))
        }
      }
      CLT[[index.CLT]] = arr.prob
    }
    names(CLT) = MAF_interval
    return(CLT)
  }
}

SPA_G.one.SNP = function(obj.precond,
                         geno,       # the order of g should be consistent with Resid and GRM!
                         Cutoff = 2,
                         impute.method = "fixed",
                         missing.cutoff = 0.15,
                         min.maf = 0.0001,
                         G.model = "Add",
                         zeta = 0,
                         tol = 1e-5)
{
  Resid = obj.precond$Resid;
  Resid.unrelated.outliers = obj.precond$Resid.unrelated.outliers;
  R_GRM_R = obj.precond$R_GRM_R;
  R_GRM_R_TwoSubjOutlier = obj.precond$R_GRM_R_TwoSubjOutlier;
  sum_R_nonOutlier = obj.precond$sum_R_nonOutlier;
  R_GRM_R_nonOutlier = obj.precond$R_GRM_R_nonOutlier;
  TwoSubj_list = obj.precond$TwoSubj_list;
  ThreeSubj_list = obj.precond$ThreeSubj_list;
  MAF_interval = obj.precond$MAF_interval;
  
  ## calculate MAF and update genotype vector
  MAF = mean(geno, na.rm = T)/2
  N = length(geno)
  pos.na = which(is.na(geno))
  missing.rate = length(pos.na)/N
  
  if(missing.rate != 0){
    if(impute.method == "fixed")
      geno[pos.na] = 2*MAF
    else
      stop("impute.method should be 'fixed'.")
  }
  
  if(MAF > 0.5){
    MAF = 1 - MAF
    geno = 2 - geno
  }
  
  if(MAF < min.maf | missing.rate > missing.cutoff)
    return(c(MAF, missing.rate, NA, NA, NA))
  
  if(G.model=="Add"){}   # do nothing if G.Model is "Add"
  if(G.model=="Dom") g = ifelse(g>=1,1,0)
  if(G.model=="Rec") g = ifelse(g<=1,0,1)
  
  ## Score statistic
  Score = sum((geno - 2*MAF) * Resid, na.rm = T)
  
  g.var = 2 * MAF * (1 - MAF)
  S.var.GRM = R_GRM_R * g.var
  z.GRM = Score/sqrt(S.var.GRM)
  
  pval.norm_GRM = pnorm(abs(z.GRM), lower.tail = FALSE)*2
  if(abs(z.GRM) < Cutoff){
    pval.spa.G_GRM = pval.norm_GRM
  }else{
    # MAF_interval = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
    order_interval = cut(MAF, breaks = MAF_interval, labels = 1:(length(MAF_interval)-1))
    order_interval = as.numeric(order_interval)
    MAF_ratio = (MAF_interval[order_interval + 1] - MAF)/(MAF_interval[order_interval + 1] - MAF_interval[order_interval])
    
    output = GetTwoProb_cpp(order_interval, MAF_ratio, R_GRM_R_TwoSubjOutlier, g.var, S.var.GRM,
                            sum_R_nonOutlier, R_GRM_R_nonOutlier, Resid.unrelated.outliers,
                            TwoSubj_list, ThreeSubj_list, Score, MAF, zeta, tol)
    
    pval.spa.G_GRM = output$pval1 + output$pval2
  }
  pval = c(pval.spa.G_GRM, pval.norm_GRM) # 2 elements
  
  output = c(MAF, missing.rate, Score, pval)
  return(output)
}
