import { QuartzConfig } from "./quartz/cfg"
import * as Plugin from "./quartz/plugins"

const config: QuartzConfig = {
  configuration: {
    pageTitle: "🛤️ Prior Paths",      
    enableSPA: true,
    enablePopovers: true,
    analytics: {
      provider: "plausible",
    },
    baseUrl: "DHintz137.github.io/quartz",
    ignorePatterns: ["private", "templates", ".obsidian"],
    defaultDateType: "created",
    theme: {
      typography: {
        header: "Schibsted Grotesk",
        body: "Source Sans Pro",
        code: "IBM Plex Mono",
      },
      colors: {
        lightMode: {             // <Key>                                           //  <defualt>
          light: "#faf8f8",      // Primary Background Color                        //  "#faf8f8"
          lightgray: "#e5e5e5",  // Secondary Background Color (search bar color)   //  "#e5e5e5"
          gray: "#b8b8b8",       // Accent Color 1  (subitle text color)            //  "#b8b8b8"
          darkgray: "#4e4e4e",   // Accent Color 2 (main text color)                //  "#4e4e4e"
          dark: "#284b63",       // Text Color (title text color)                   //  "#284b63"
          secondary: "#284b63",  // Highlight/Detail Color (link text color)        //  "#284b63"
          tertiary: "#84a59d",
          highlight: "rgba(143, 159, 169, 0.15)",
        },
        darkMode: {             // <Key>                                             // <defualt>
          light: "#2E2E2E",     // Primary Background Color                          //  defualt: #161618
          lightgray: "#393639", // Secondary Background Color (search bar color)     //  defualt: #393639
          gray: "#F1B923",      // Accent Color 1  (subitle text, read time color)  
          darkgray: "#d4d4d4",  // Accent Color 2 (main text color)  
          dark: "#ebebec",      // Text Color (title text color)   
          secondary: "#F6F4D6", // Highlight/Detail Color (link text/ website name color)  
          tertiary: "#84a59d",
          highlight: "rgba(143, 159, 169, 0.15)",  // or is hex its #8F9FA9
        },
      },
    },
  },
  plugins: {
    transformers: [
      Plugin.FrontMatter(),
      Plugin.TableOfContents(),
      Plugin.CreatedModifiedDate({
        priority: ["frontmatter", "filesystem"], // you can add 'git' here for last modified from Git but this makes the build slower
      }),
      Plugin.SyntaxHighlighting(),
      Plugin.ObsidianFlavoredMarkdown({ enableInHtmlEmbed: false }),
      Plugin.GitHubFlavoredMarkdown(),
      Plugin.CrawlLinks({ markdownLinkResolution: "shortest" }),
      Plugin.Latex({ renderEngine: "katex" }),
      Plugin.Description(),
    ],
    filters: [Plugin.RemoveDrafts()],
    emitters: [
      Plugin.AliasRedirects(),
      Plugin.ComponentResources({ fontOrigin: "googleFonts" }),
      Plugin.ContentPage(),
      Plugin.FolderPage(),
      Plugin.TagPage(),
      Plugin.ContentIndex({
        enableSiteMap: true,
        enableRSS: true,
      }),
      Plugin.Assets(),
      Plugin.Static(),
      Plugin.NotFoundPage(),
    ],
  },
}

export default config
