using ArcGIS.Desktop.Core;
using ArcGIS.Desktop.Framework;
using ArcGIS.Desktop.Framework.Contracts;


namespace Combat_Search_and_Rescue_Tools
{
    internal class CSARToolsModule : Module
    {
        private static CSARToolsModule _this = null;

        /// <summary>
        /// Retrieve the singleton instance to this module here
        /// </summary>
        public static CSARToolsModule Current
        {
            get
            {
                return _this ?? (_this = (CSARToolsModule)FrameworkApplication.FindModule("Combat_Search_and_Rescue_Tools_Module"));
            }
        }

        #region Overrides
        /// <summary>
        /// Called by Framework when ArcGIS Pro is closing
        /// </summary>
        /// <returns>False to prevent Pro from closing, otherwise True</returns>
        protected override bool CanUnload()
        {
            //TODO - add your business logic
            //return false to ~cancel~ Application close
            return true;
        }

        #endregion Overrides
    }
}
