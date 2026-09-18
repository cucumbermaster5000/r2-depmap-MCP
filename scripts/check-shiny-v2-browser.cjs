// Set NODE_PATH to a local Playwright installation; no package download needed.
const { chromium } = require('playwright');
const fs = require('node:fs');
const path = require('node:path');
(async () => {
  const folder = path.resolve('outputs/shiny-v2-acceptance/browser');
  fs.mkdirSync(folder, {recursive:true});
  const browser = await chromium.launch({channel:'msedge',headless:true});
  try {
    const page = await browser.newPage({viewport:{width:1550,height:1100},acceptDownloads:true});
    const errors=[]; page.on('pageerror', e=>errors.push(e.message));
    await page.goto(process.env.SHINY_TEST_URL || 'http://127.0.0.1:3878');
    await page.waitForFunction(()=>window.Shiny && Shiny.shinyapp && Shiny.shinyapp.$socket.readyState===1);
    for(const gene of ['HLX','MYC','SLC2A1']) {
      await page.locator('#gene').fill(gene);
      await page.locator('#run').click();
      await page.waitForFunction(g=>document.querySelector('#status')?.innerText.includes('Analysis complete: '+g),gene,{timeout:180000});
      await page.getByRole('tab',{name:'MB Subtypes',exact:true}).click();
      await page.waitForFunction(()=>document.querySelectorAll('#subtype_pairs tbody tr').length===66);
      await page.waitForFunction(()=>document.querySelector('#subtype_plot img')?.naturalWidth>0);
      const summary=await page.locator('#subtype_summary').innerText();
      if(!summary.includes('Annotated: 763') || !summary.includes('Subtypes: 12') || !summary.includes('Missing subtype: 0')) throw Error(summary);
      if(await page.locator('#subtype_descriptives tbody tr').count()!==12) throw Error('Wrong descriptive row count');
      if(!(await page.locator('#subtype_overall').innerText()).includes('Kruskal')) throw Error('Missing overall test');
      await page.screenshot({path:path.join(folder,gene+'-subtypes.png'),fullPage:true});
      await page.getByRole('tab',{name:'MB Subgroups',exact:true}).click();
      await page.waitForFunction(()=>document.querySelector('#subgroup_plot img')?.naturalWidth>0);
      await page.getByRole('tab',{name:'Downloads',exact:true}).click();
      for(const id of ['patients','stats_csv','figure_pdf','figure_png','subtype_patients','subtype_summaries','subtype_pairwise','subtype_overall_csv','subtype_pdf','subtype_png']) {
        const pending=page.waitForEvent('download');
        await page.locator('#'+id).click();
        const download=await pending;
        if(await download.failure()) throw Error(await download.failure());
        const destination=path.join(folder,download.suggestedFilename());
        await download.saveAs(destination);
        if(fs.statSync(destination).size<100) throw Error('Empty download '+id);
      }
      console.log(gene+': 12 subtypes, 66 comparisons, both plots, 10 browser downloads passed');
    }
    if(errors.length) throw Error(errors.join('\n'));
    fs.writeFileSync(path.join(folder,'acceptance.json'),JSON.stringify({genes:['HLX','MYC','SLC2A1'],downloads:30,javascriptErrors:errors},null,2));
  } finally { await browser.close(); }
})().catch(e=>{console.error(e);process.exitCode=1;});
