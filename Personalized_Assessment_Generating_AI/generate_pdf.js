const puppeteer = require('puppeteer');
const fs = require('fs');

(async () => {
    const browser = await puppeteer.launch();
    const page = await browser.newPage();

    // Load the HTML file with MathJax included
    const htmlPath = process.argv[2];
    const outputPdf = process.argv[3];
    const content = fs.readFileSync(htmlPath, 'utf-8');

    await page.setContent(content, { waitUntil: 'networkidle0' });

    // Generate the PDF
    await page.pdf({ path: outputPdf, format: 'A4' });

    await browser.close();
})();