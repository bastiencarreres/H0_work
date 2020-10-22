def chi2(v):
    n = len(sn_host)
    mu_mu4258=v[:n]
    zp_w4258=v[n]
    b_w = v[n+1]
    m_v4258 = v[n+2]
    Z_w = v[n+3]
    chi2 = 0
    for host,i in zip(sn_host,np.arange(n)):
        selec_sn = data_sn[data_sn['Field']==host]
        m_sn_mod = m_v4258 + mu_mu4258[i]
        m_sn = selec_sn['mag_5a'] - 5*av
        cov_sn = pw(selec_sn['sigma_mag_5a'],2)-25*pw(sigma_av,2)
        chi2 += pw(m_sn-m_sn_mod,2)/cov_sn
        
        selec_ceph = data_cephe[data_cephe['Field']== host]
        C_ceph_inv=np.linalg.inv(np.diag(pw(selec_ceph['sigma_tot'],2)))
        m_ceph_mod = mu_mu4258[i] + zp_w4258 + b_w*np.log10(selec_ceph['Per']) + Z_w*(selec_ceph['[O/H]']-12)
        m_ceph = selec_ceph['F160W']-R*selec_ceph['F555W-F814W']
        chi2 += (m_ceph-m_ceph_mod).T @ C_ceph_inv @ (m_ceph-m_ceph_mod)
        
    ceph_4258 = data_cephe[data_cephe['Field']==b'N4258']
    m_ceph_mod = zp_w4258 + b_w*np.log10(ceph_4258['Per']) + Z_w*(ceph_4258['[O/H]']-12)
    m_ceph = ceph_4258['F160W']-R*ceph_4258['F555W-F814W']
    C_ceph_inv=np.linalg.inv(np.diag(pw(ceph_4258['sigma_tot'],2)))
    chi2 += (m_ceph-m_ceph_mod).T @ C_ceph_inv @ (m_ceph-m_ceph_mod)
    
    return chi2

