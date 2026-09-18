const {chromium}=require('playwright');
const fs=require('node:fs');const path=require('node:path');
(async()=>{
 const folder=path.resolve('outputs/shiny-v3-acceptance/browser');fs.mkdirSync(folder,{recursive:true});
 const browser=await chromium.launch({channel:'msedge',headless:true});
 try {
  const page=await browser.newPage({viewport:{width:1550,height:1050},acceptDownloads:true});
  const errors=[];page.on('pageerror',e=>errors.push(e.message));let downloads=0;
  await page.goto(process.env.SHINY_TEST_URL||'http://127.0.0.1:3879');
  await page.waitForFunction(()=>window.Shiny?.shinyapp?.$socket?.readyState===1);
  for(const gene of ['HLX','MYC','SLC2A1']) {
   await page.locator('#gene').fill(gene);await page.locator('#run').click();
   await page.waitForFunction(g=>document.querySelector('#status')?.innerText.includes('Analysis complete: '+g),gene,{timeout:180000});
   await page.getByRole('tab',{name:'Overview',exact:true}).click();
   await page.waitForFunction(()=>document.querySelector('#overview')?.innerText.includes('Clinical results'));
   const overview=await page.locator('#overview').innerText();
   if(!overview.includes('survival N=612') || !overview.includes('metastatic N=176')) throw Error('Missing overview results');
   const cohorts=gene==='MYC'?['all']:['all','group3'];
   for(const ep of ['metastasis','survival']) {
    await page.getByRole('tab',{name:ep[0].toUpperCase()+ep.slice(1),exact:true}).click();
    await page.waitForFunction(({ep,gene})=>document.querySelector('#'+ep+'_content h3')?.innerText===gene+' '+ep,{ep,gene});
    for(const co of cohorts) {
     const id=ep+'_'+co;
     await page.waitForFunction(id=>document.querySelector('#'+id+'_plot img')?.naturalWidth>0,id);
     if((await page.locator('#'+id+'_stats').innerText()).length<100)throw Error('Missing statistics '+id);
    }
    if(gene==='MYC' && await page.locator('#'+ep+'_group3_plot').count())throw Error('Unselected subgroup displayed');
    await page.locator('#'+ep+'_all_plot').screenshot({path:path.join(folder,gene+'-'+ep+'.png')});
   }
   for(const [tab,id] of [['MB Subgroups','subgroup_plot'],['MB Subtypes','subtype_plot']]) {
    await page.getByRole('tab',{name:tab,exact:true}).click();
    await page.waitForFunction(id=>document.querySelector('#'+id+' img')?.naturalWidth>0,id);
   }
   await page.getByRole('tab',{name:'Downloads',exact:true}).click();
   const ids=['patients','stats_csv','figure_pdf','figure_png','subtype_patients','subtype_summaries','subtype_pairwise','subtype_overall_csv','subtype_pdf','subtype_png'];
   for(const ep of ['metastasis','survival'])for(const co of cohorts)for(const kind of ['patient_data','statistics','plot_pdf','plot_png'])ids.push(ep+'_'+co+'_download_'+kind);
   for(const id of ids) {
    const pending=page.waitForEvent('download');await page.locator('#'+id).click();
    const download=await pending;if(await download.failure())throw Error(await download.failure());
    const destination=path.join(folder,download.suggestedFilename());await download.saveAs(destination);
    if(fs.statSync(destination).size<50)throw Error('Empty download');downloads++;
   }
   console.log(gene+': Overview, V1/V2 plots, clinical plots, selected cohorts and '+ids.length+' downloads passed');
  }
  if(errors.length)throw Error(errors.join('\n'));
  fs.writeFileSync(path.join(folder,'acceptance.json'),JSON.stringify({genes:['HLX','MYC','SLC2A1'],downloads,javascriptErrors:errors},null,2));
 }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
